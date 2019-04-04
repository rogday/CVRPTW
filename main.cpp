#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <list>
#include <unordered_set>

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

  struct path_t {
    std::vector<u64> customers;
    i64 capacity = 0;
    i64 time = 0;
  };

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
      std::vector<sf::Vertex> line(path.customers.size());

      int r = prng() % 256, g = prng() % 256, b = prng() % 256;
      sf::Color color(r, g, b);

      int i = 0;
      for (auto &customer : path.customers) {
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

      i64 current_time = 0;
      i64 capacity = vehicle_capacity;

      path.customers.push_back(0);

      auto insert = [&path, &current_time, &capacity, &set,
                     this](i64 next, bool debug = false) {
        current_time += i64(distance(path.customers.back(), next) + 0.5);
        current_time +=
            std::max(customers[next].ready_time - current_time, 0ll);

        if (next != 0) {
          if (capacity >= customers[next].demand)
            set.erase(next);
          else
            customers[next].demand -= capacity;

          capacity -= std::min(capacity, customers[next].demand);
        } else
          capacity = vehicle_capacity;

        std::cout << next << " ";

        if (debug)
          std::cout << ", time: " << current_time
                    << ", distance: " << distance(path.customers.back(), next)
                    << ", capacity: " << capacity
                    << ", i=" << path.customers.back() << ", j=" << next
                    << std::endl;

        path.customers.push_back(next);
      };

      while (current_time < warehouse_closing) {
        i64 next = -1;
        double max_profit = 0.0;

        for (u64 i : set) {
          i64 open_time = customers[i].ready_time - current_time;

          if (distance(path.customers.back(), i) >
              warehouse_closing + distance(path.customers.back(), 0)) {
            next = -1;
            break;
          }

          double profit = profitability(path.customers.back(), i, open_time,
                                        customers[i].due_date - current_time -
                                            customers[i].service_time,
                                        capacity);

          // std::cout << profit << std::endl;
          if (profit > max_profit) {
            max_profit = profit;
            next = i;
          }
        }

        if (next != -1 && current_time < customers[0].due_date) {
          insert(next);
        } else
          break;
      }

      if (path.customers.back() != 0)
        insert(0);

      path.time = current_time;
      path.capacity = capacity;

      std::cout << std::endl << std::endl;

      if (set.size() == 1)
        break;
    }
  }

  void perturbation(u64 path) {
    double acc = 0;
    for (u64 i = 0, k = 0; i < paths.size(); ++i) {
      if (i != path && !paths[i].customers.empty()) {
        while (k != paths[i].customers.size() - 1) {
          double t = estimate_time(paths[i], k);
          acc += t;
          std::cout << t << std::endl;
        }
        k = 0;
      }
    }
    std::cout << acc << std::endl;
  }

  double estimate_time(path_t &path, u64 &start) {
    double time = 0.0;

    std::cout << std::endl << std::endl;

    do {
      std::cout << path.customers[start] << " ";
      time += distance(path.customers[start], path.customers[start + 1]);
      ++start;
    } while (path.customers[start] != 0);

    std::cout << path.customers[start] << std::endl << std::endl;
    time += distance(path.customers[start], 0);

    return time;
  }

  double distance(u64 i, u64 j) {
    return std::hypot(customers[i].x - customers[j].x,
                      customers[i].y - customers[j].y);
  }

  // ->> max
  double profitability(u64 start, u64 end, i64 open_time, i64 close_time,
                       i64 capacity) {
    double d = distance(start, end);
    double ret = 0.0;
    open_time = (open_time <= 0) ? 1 : open_time;

    if (d == 0)
      return 0;

    if (capacity >= vehicle_capacity / 3 && end == 0)
      return 0;

    if (close_time - d > 0) {
      ret = 1.0 / (open_time) + 1.0 * max_distance / d;

      if (capacity > customers[end].demand)
        ret += 1.0;
      else {
        if (customers[end].demand != 0)
          ret += 1.0 * capacity / customers[end].demand;
        else if (capacity == 0)
          return 100;
      }
    }

    return ret;
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
