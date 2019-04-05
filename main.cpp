#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include <filesystem>

#include <SFML/Graphics.hpp>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using path_map_t = std::vector<std::string>;

class vrp_data_storage {
  struct customer_description_t {

    i64 identificator, x, y, demand, ready_time, due_date, service_time;

    customer_description_t(i64 id, i64 x, i64 y, i64 d, i64 rt, i64 dd, i64 st)
        : identificator{id}, x{x}, y{y}, demand{d},
          ready_time{rt}, due_date{dd}, service_time{st} {}

    friend std::ostream &operator<<(std::ostream &o,
                                    customer_description_t const &customer) {
      return o << "customer{\nid:\t\t" << customer.identificator << "\n"
               << "x: \t\t" << customer.x << "\n"
               << "y: \t\t" << customer.y << "\n"
               << "demand: \t" << customer.demand << "\n"
               << "ready time: \t" << customer.ready_time << "\n"
               << "service time: \t" << customer.due_date << "\n}\n";
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
    find_world_bounds();

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

  void local_search() {}

  double overall_distance() {
    double ret = 0.0;

    for (auto &path : paths)
      ret += estimate_time(path);

    return ret;
  }

private:
  void clear() {
    max_x = max_y = std::numeric_limits<i64>::min();
    min_x = min_y = std::numeric_limits<i64>::max();

    vehicle_amount = vehicle_capacity = 0;
    customers.clear();
    paths.clear();
  }

  std::pair<i64, i64> find_best_neighbour(std::unordered_set<u64> &set,
                                          path_t &path, i64 current_time,
                                          i64 capacity) {
    i64 next = -1, warehouse_closing = customers[0].due_date;
    double min_time = std::numeric_limits<double>::max();

    for (u64 i : set) {
      // distance from current to new + current time + await time
      double time = current_time + distance(path.back(), i);
      time += std::max(customers[i].ready_time - time, 0.0);

      double time_served = time + customers[i].service_time;

      // if we CAN go to the new location and not end up getting late for
      // warehouse closing and we can fullfill the needs of a client
      if (time_served <= warehouse_closing - distance(i, 0) &&
          time <= customers[i].due_date && capacity >= customers[i].demand &&
          time_served < min_time) {
        min_time = time_served;
        next = i;
      }
    }

    return {next, min_time};
  }

  void greedy_heuristics() {
    std::unordered_set<u64> set;
    for (u64 i = 0; i < customers.size(); ++i)
      set.insert(i);

    i64 warehouse_closing = customers[0].due_date;

    for (auto &path : paths) {
      i64 current_time = 0;
      i64 capacity = vehicle_capacity;

      path.push_back(0);
      while (current_time < warehouse_closing) {
        auto [next, min_time] =
            find_best_neighbour(set, path, current_time, capacity);

        if (next != -1) { // insertion
          current_time = min_time;
          set.erase(next);
          capacity -= customers[next].demand;
          path.push_back(next);
        } else
          break;
      }

      if (path.size() != 1) // heading back to the depot
        path.push_back(0);
      else
        break;
    }

    assert(set.empty());
  }

  void perturbation(u64 path) {}

  double estimate_time(path_t &path) {
    double time = 0.0;

    for (std::size_t i = 0; i + 1 < path.size(); ++i)
      time += distance(path[i], path[i + 1]);

    return time;
  }

  double distance(u64 i, u64 j) {
    return std::hypot(customers[i].x - customers[j].x,
                      customers[i].y - customers[j].y);
  }

  void find_world_bounds() {
    std::uint64_t size = customers.size();

    for (std::size_t i = 0; i < size; ++i)
      for (std::size_t j = i; j < size; ++j) {

        max_x = std::max(max_x, customers[i].x);
        max_y = std::max(max_y, customers[i].y);

        min_x = std::min(min_x, customers[i].x);
        min_y = std::min(min_y, customers[i].y);
      }
  }

  u64 vehicle_amount;
  u64 vehicle_capacity;
  std::vector<customer_description_t> customers;

  i64 max_x, min_x;
  i64 max_y, min_y;

  std::vector<path_t> paths;
};

sf::RenderWindow &init_window(double size) {
  sf::VideoMode vm = sf::VideoMode::getDesktopMode();

  static sf::RenderWindow window(
      sf::VideoMode(int(vm.width / size), int(vm.height / size)), "CVRPTW",
      sf::Style::Titlebar | sf::Style::Close | sf::Style::Resize);

  window.setPosition(sf::Vector2i(vm.width / 2 - vm.width / (size * 2),
                                  vm.height / 2 - vm.height / (size * 2)));

  window.setKeyRepeatEnabled(true);
  window.setVerticalSyncEnabled(true);

  return window;
}

auto get_map(std::string path) {
  path_map_t map;

  for (auto entry : std::filesystem::directory_iterator(path))
    map.push_back(entry.path().string());

  return map;
}

void print_choice(path_map_t &map) {
  std::size_t i = -1;
  std::cout << std::endl;
  while (++i != map.size())
    std::cout << "#" << i << ": " << map[i] << std::endl;
  std::cout << std::endl;
}

int main() {
  vrp_data_storage vrp;

  vrp.read_data("..\\input\\C108.txt");
  vrp.generate_initial_solution(vrp_data_storage::Heuristics::Greedy);
  vrp.local_search();

  auto &window = init_window(1.5);

  auto regular_map = get_map("..\\input");
  auto bonus_map = get_map("..\\bonus");
  std::reference_wrapper<path_map_t> current_map = regular_map;

  bool bonus = false;
  sf::Event event;
  while (window.isOpen()) {
    window.clear(sf::Color::Black);
    vrp.draw(window);

    while (window.pollEvent(event)) {
      switch (event.type) {

      case sf::Event::KeyPressed:
        if (event.key.code == sf::Keyboard::Escape) { // change colors
          window.close();
        }
        if (event.key.code == sf::Keyboard::Space) { // change colors
          vrp.draw(window, true);

        } else if (event.key.code == sf::Keyboard::P) { // print
          print_choice(current_map);
        } else if (event.key.code == sf::Keyboard::B) { // flip bonus

          bonus ^= true;
          current_map = bonus ? bonus_map : regular_map;
          std::cout << (bonus ? "bonus tasks" : "regular tasks") << std::endl;

        } else if (event.key.code >= sf::Keyboard::Num0 &&
                   event.key.code <= sf::Keyboard::Num9) {

          std::size_t n = event.key.code - sf::Keyboard::Num0;
          if (n >= current_map.get().size())
            break;

          std::string path = current_map.get()[n];
          vrp.read_data(path);
          vrp.generate_initial_solution(vrp_data_storage::Heuristics::Greedy);

          std::cout << "#" << n << " "
                    << std::filesystem::path(path).stem().string() << ": "
                    << vrp.overall_distance() << std::endl;
        }

        break;

      case sf::Event::Closed:
        window.close();
        break;

      case sf::Event::Resized:
        sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
        window.setView(sf::View(visibleArea));
        break;
      }
    }
  }

  return 0;
}