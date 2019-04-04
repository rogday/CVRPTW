#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <list>
#include <unordered_set>
#include <cassert>

#include <SFML/Graphics.hpp>

using u64 = std::uint64_t;
using i64 = std::int64_t;

namespace {

template <typename T, typename U>
std::ostream &operator<<(std::ostream &o, std::pair<T, U> const &pair) {
  return o << pair.first << " " << pair.second;
}

template <typename T>
std::ostream &operator<<(std::ostream &o, std::vector<T> const &vector) {
  for (auto const &element : vector)
    o << element << std::endl << std::endl;
  return o;
}

} // namespace

class vrp_data_storage {
  struct customer_description_t {

    i64 identificator, x, y, demand, ready_time, due_date, service_time;

    customer_description_t(i64 id, i64 x, i64 y, i64 d, i64 rt, i64 dd, i64 st)
        : identificator{id}, x{x}, y{y}, demand{d},
          ready_time{rt}, due_date{dd}, service_time{st} {}

    friend std::ostream &operator<<(std::ostream &o,
                                    customer_description_t const &customer) {
      return o << "NO:\t\t" << customer.identificator << "\n"
               << "XC: \t\t" << customer.x << "\n"
               << "YC: \t\t" << customer.y << "\n"
               << "DEMAND: \t" << customer.demand << "\n"
               << "READY TIME: \t" << customer.ready_time << "\n"
               << "SERVICE TIME: \t" << customer.due_date << "\n";
    }
  };

  using path_t = std::vector<u64>;

public:
  enum Heuristics { Greedy = 0 };

  vrp_data_storage() = default;

  void read_data(std::string const &filename) {
    this->clear();

    std::fstream fstream{filename, std::ios_base::in};
    std::string current_line;

    while (std::getline(fstream, current_line)) {
      std::istringstream iss(current_line);
      if (iss >> vehicle_amount >> vehicle_capacity) {
        paths.resize(vehicle_amount);
        break;
      }
    }

    i64 id, x, y, demand, ready_time, due_date, service_time;
    while (std::getline(fstream, current_line)) {
      std::istringstream iss(current_line);
      if (iss >> id >> x >> y >> demand >> ready_time >> due_date >>
          service_time)
        customers.emplace_back(id, x, y, demand, ready_time, due_date,
                               service_time);
    }
  }

  void generate_initial_solution(Heuristics heuristics = Heuristics::Greedy) {
    generate_distance_matrix();

    switch (heuristics) {
    case Heuristics::Greedy:
      greedy_heuristics();
      break;
    default:
      std::cerr << "Wrong type." << std::endl;
      break;
    }
  }

  void draw(sf::RenderWindow &window, bool new_colors = false) {
    static std::random_device rd;
    static unsigned int seed = rd();

    if (new_colors)
      seed = rd();

    std::mt19937 prng(seed);

    sf::Vector2u window_size = window.getSize();
    double zoom = 0.9;

    auto world_to_window = [window_size, zoom,
                            this](i64 current_x,
                                  i64 current_y) -> sf::Vector2<i64> {
      i64 x =
          (current_x - min_x) / double(max_x - min_x) * window_size.x * zoom +
          (1.0 - zoom) * window_size.x / 2;
      i64 y =
          (current_y - min_y) / double(max_y - min_y) * window_size.y * zoom +
          (1.0 - zoom) * window_size.y / 2;
      return {x, y};
    };

    for (auto &path : paths) {
      std::vector<sf::Vertex> line(path.size());

      int r = prng() % 256, g = prng() % 256, b = prng() % 256;
      sf::Color color(r, g, b);

      int i = 0;
      for (auto &customer : path) {
        auto [x, y] =
            world_to_window(customers[customer].x, customers[customer].y);

        line[i].position = sf::Vector2f(x, y);
        line[i].color = color;
        ++i;
      }
      window.draw(line.data(), line.size(), sf::LineStrip);
    }
    window.display();
  }

  void local_search() { perturbation(12); }

private:
  void clear() {
    max_x = max_y = std::numeric_limits<i64>::min();
    min_x = min_y = std::numeric_limits<i64>::max();

    vehicle_amount = vehicle_capacity = 0;
    customers.clear();
    max_distance = 0.0;
    distance_matrix.clear();
    paths.clear();
  }

  void greedy_heuristics() {
    std::unordered_set<u64> set;
    for (u64 i = 0; i < customers.size(); ++i)
      set.insert(i);

    i64 warehouse_closing = customers[0].due_date;

    for (auto &path : paths) {

      i64 current_time = 0, new_time = 0;
      i64 capacity = vehicle_capacity;

      path.push_back(0);

      while (current_time < warehouse_closing) {
        i64 next = -1;
        double min_dist = std::numeric_limits<double>::max();

        for (u64 i : set) {
          // distance from current to new + current time + time of service +
          // await time
          i64 tmp = current_time + distance(path.back(), i);
          tmp += std::max(customers[i].ready_time - tmp, 0ll) +
                 customers[i].service_time;

          // if we CAN go to the new location and not end up getting late for
          // warehouse closing and we can fullfill the needs of a client
          if (tmp <= warehouse_closing - distance(i, 0) &&
              tmp <= customers[i].due_date && capacity >= customers[i].demand) {
            double dist = distance(path.back(), i);

            if (dist < min_dist) {
              new_time = tmp;
              min_dist = dist;
              next = i;
              break;
            }
          }
        }

        if (next != -1) { // insertion
          current_time = new_time;
          set.erase(next);
          capacity -= customers[next].demand;
          path.push_back(next);

          if (false) {
            std::cout << ", time: " << current_time
                      << ", distance: " << distance(path.back(), next)
                      << ", capacity: " << capacity << ", i=" << path.back()
                      << ", j=" << next << std::endl;
          }
        } else
          break;
      }

      if (path.size() > 1) // heading back to the depot
        path.push_back(0);

      if (set.empty())
        break;
    }
    assert(set.empty());
  }

  void perturbation(u64 path) {
    double acc = 0;
    int s = 1;

    for (u64 i = 0, k = 0; i < paths.size(); ++i) {
      if (i != path && !paths[i].empty()) {
        ++s;
        while (k != paths[i].size() - 1) {
          double t = estimate_time(paths[i], k);
          acc += t;
        }
        k = 0;
      }
    }
    std::cout << acc << " vehicles: " << s << std::endl;
  }

  double estimate_time(path_t &path, u64 &start) {
    double time = 0.0;

    do {
      time += distance(path[start], path[start + 1]);
      ++start;
    } while (path[start] != 0);

    time += distance(path[start], 0);

    return time;
  }

  double distance(u64 i, u64 j) {
    return std::hypot(customers[i].x - customers[j].x,
                      customers[i].y - customers[j].y);
  }

  std::pair<u64, u64> generate_distance_matrix() {
    std::uint64_t size = customers.size();

    distance_matrix.resize(size);
    std::fill(distance_matrix.begin(), distance_matrix.end(),
              std::vector<double>(size, 0.0));

    u64 start = 0, end = 0;

    for (std::size_t i = 0; i < size; ++i)
      for (std::size_t j = i; j < size; ++j) {
        double d = std::sqrt(std::pow(customers[i].x - customers[j].x, 2) +
                             std::pow(customers[i].y - customers[j].y, 2));
        distance_matrix[i][j] = distance_matrix[j][i] = d;

        max_x = std::max(max_x, customers[i].x);
        max_y = std::max(max_y, customers[i].y);

        min_x = std::min(min_x, customers[i].x);
        min_y = std::min(min_y, customers[i].y);

        if (d > max_distance) {
          max_distance = d;
          start = i;
          end = j;
        }
      }

    return {start, end};
  }

  u64 vehicle_amount;
  u64 vehicle_capacity;
  std::vector<customer_description_t> customers;

  double max_distance;
  i64 max_x, min_x;
  i64 max_y, min_y;

  std::vector<std::vector<double>> distance_matrix;

  std::vector<path_t> paths;
};

int main() {
  srand(time(nullptr));

  vrp_data_storage vrp;
  vrp.read_data("..\\input\\C108.txt");
  vrp.generate_initial_solution(vrp_data_storage::Heuristics::Greedy);

  vrp.local_search();

  using namespace sf;

  VideoMode vm = VideoMode::getDesktopMode();

  double size = 1.5;

  RenderWindow window(VideoMode(int(vm.width / size), int(vm.height / size)),
                      "Test.", Style::Titlebar | Style::Close | Style::Resize);

  window.setPosition(Vector2i(vm.width / 2 - vm.width / (size * 2),
                              vm.height / 2 - vm.height / (size * 2)));

  window.setKeyRepeatEnabled(true);
  window.setVerticalSyncEnabled(true);

  Event event;
  while (window.isOpen()) {
    vrp.draw(window);

    while (window.pollEvent(event)) {
      switch (event.type) {

      case Event::KeyPressed:
        if (event.key.code == Keyboard::Space)
          vrp.draw(window, true);

        break;

      case Event::Closed:
        window.close();
        break;

      case Event::Resized:
        FloatRect visibleArea(0, 0, event.size.width, event.size.height);
        window.setView(View(visibleArea));
        break;
      }
    }
  }

  return 0;
}