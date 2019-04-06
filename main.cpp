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
  using iter_t = path_t::iterator;

public:
  enum Heuristics { Greedy = 0 };

  vrp_data_storage() = default;
  vrp_data_storage(sf::RenderWindow &window) : window(window) {}

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

    sanity_check();
  }

  void draw(bool new_colors = false) {
    static std::random_device rd;
    static unsigned int seed = rd();

    window.clear(sf::Color::Black);

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

  void perform_local_search() { local_search(); }

  double overall_distance() {
    double ret = 0.0;
    for (auto &path : paths)
      ret += estimate_time(path);
    return ret;
  }

  void sanity_check() {
    fill_set();
    bool is_ok = true;
    for (auto &path : paths)
      is_ok &= valid_path(path);
    if (!is_ok)
      std::cerr << "sanity_check has failed." << std::endl;
  }

private:
  void clear() {
    max_x = max_y = std::numeric_limits<i64>::min();
    min_x = min_y = std::numeric_limits<i64>::max();

    vehicle_amount = vehicle_capacity = 0;
    customers.clear();
    paths.clear();
  }

  void fill_set() {
    set.clear();
    for (u64 i = 1; i < customers.size(); ++i)
      set.insert(i);
  }

  double can_be_neighbour(double current_time, u64 old_one, u64 new_one,
                          i64 capacity) {
    i64 warehouse_closing = customers[0].due_date;

    // distance from current to new + current time + await time
    double time = current_time + dist[old_one][new_one];
    time += std::max(customers[new_one].ready_time - time, 0.0);

    double time_served = time + customers[new_one].service_time;

    // if we CAN go to the new location and not end up getting late for
    // warehouse closing and we can fullfill the needs of a client
    if (time_served <= warehouse_closing - dist[new_one][0] &&
        time <= customers[new_one].due_date &&
        capacity >= customers[new_one].demand)
      return time_served;

    return -1.0;
  }

  bool update_variables(u64 index, double &current_time, double new_time,
                        i64 &capacity) {
    current_time = new_time;
    bool ret = ((set.erase(index) == 1) || (index == 0));
    capacity -= customers[index].demand;
    return ret;
  }

  void insert_customer_to_path(path_t &path, u64 index, double &current_time,
                               double new_time, i64 &capacity) {
    update_variables(index, current_time, new_time, capacity);
    path.push_back(index);
  }

  std::pair<i64, double> find_best_neighbour(path_t &path, double current_time,
                                             i64 capacity) {
    i64 next = -1;
    double min_time = std::numeric_limits<double>::max();

    for (u64 i : set) {
      double time = can_be_neighbour(current_time, path.back(), i, capacity);

      if (time >= 0 && time < min_time) {
        min_time = time;
        next = i;
      }
    }

    return {next, min_time};
  }

  void greedy_heuristics() {
    fill_set();

    u64 count = 0;
    for (auto &path : paths) {
      double current_time = 0;
      i64 capacity = vehicle_capacity;

      path.push_back(0);
      while (true) {
        auto [next, min_time] =
            find_best_neighbour(path, current_time, capacity);

        if (next != -1)
          insert_customer_to_path(path, next, current_time, min_time, capacity);
        else
          break;
      }

      if (path.size() != 1) { // heading back to the depot
        path.push_back(0);
        ++count;
      } else {
        paths.resize(count);
        break;
      }
    }

    if (!set.empty())
      std::cerr << "customers set is not empty" << std::endl;
  }

  bool valid_path(path_t &path) {
    double current_time = 0;
    i64 capacity = vehicle_capacity;

    for (auto it = std::begin(path); it != std::prev(std::end(path)); ++it) {
      double time =
          can_be_neighbour(current_time, *it, *std::next(it), capacity);
      if (time < 0)
        return false;
      if (!update_variables(*std::next(it), current_time, time, capacity))
        return false;
    }

    return true;
  }

  // a -> b -> c   >  a -> c -> b
  // everyone says it's garbage, so I'll skip it for now
  // void two_opt(){};

  // a -> b -> c   >  a -> b -> e -> c
  // d -> e -> f   >  d -> f
  void relocate(path_t &p1, path_t &p2, iter_t to, iter_t &from) {
    p1.insert(to, *from);
    from = p2.erase(from);
  }

  // a -> b -> c   >  a -> e -> c
  // d -> e -> f   >  d -> f -> b
  void exchange(iter_t to, iter_t from) { std::swap(*to, *from); }

  // a -> b -> c   >  a -> e -> f
  // d -> e -> f   >  d -> b -> c
  void cross(path_t &p1, path_t &p2, iter_t in, iter_t from, iter_t to) {
    // p1.splice(in, p2, from, to);
  }

  void local_search() {
    fill_set();
    double current_distance, min_distance = overall_distance(), last_distance;

    do {
      last_distance = min_distance;

      for (u64 i = 0; i < std::size(paths); ++i)
        for (u64 j = i + 1; j < std::size(paths); ++j) {
          auto &p1 = paths[i];
          auto &p2 = paths[j];

          path_t best1 = p1, best2 = p2;

          current_distance =
              min_distance - estimate_time(p1) - estimate_time(p2);

          for (iter_t it = std::next(std::begin(p1));
               it != std::prev(std::end(p1)); ++it)
            for (iter_t jt = std::next(std::begin(p2));
                 jt != std::prev(std::end(p2)); ++jt) {

              path_t new_p1(std::begin(p1), it);
              path_t new_p2(std::begin(p2), jt);

              new_p1.insert(std::end(new_p1), jt, std::end(p2));
              new_p2.insert(std::end(new_p2), it, std::end(p1));

              if (valid_path(new_p1) && valid_path(new_p2)) {
                double tmp = current_distance + estimate_time(new_p1) +
                             estimate_time(new_p2);

                if (tmp < min_distance) {
                  min_distance = tmp;
                  std::cout << "new_min: " << min_distance << std::endl;

                  best1 = std::move(new_p1);
                  best2 = std::move(new_p2);
                  draw();
                }
              }
              for (auto c : p1)
                set.insert(c);
              for (auto c : p2)
                set.insert(c);
            }
          p1 = best1;
          p2 = best2;
        }
      std::cout << (min_distance / last_distance) << std::endl;
    } while (1 - min_distance / last_distance > 0.0001);

    if (min_distance != overall_distance())
      std::cerr << "min_distance: " << min_distance
                << ", overall(): " << overall_distance() << std::endl;
  }

  void perturbation(u64 path) {}

  double estimate_time(path_t &path) {
    double time = 0.0;

    for (auto it = std::begin(path); it != std::prev(std::end(path)); ++it)
      time += dist[*it][*std::next(it)];

    return time;
  }

  /*
    double distance(u64 i, u64 j) {
      return std::hypot(customers[i].x - customers[j].x,
                        customers[i].y - customers[j].y);
    }*/

  void find_world_bounds() {
    std::uint64_t size = customers.size();
    dist.resize(size);
    std::fill(std::begin(dist), std::end(dist), std::vector<double>(size, 0.0));

    for (std::size_t i = 0; i < size; ++i)
      for (std::size_t j = i; j < size; ++j) {
        dist[i][j] = dist[j][i] = std::hypot(customers[i].x - customers[j].x,
                                             customers[i].y - customers[j].y);

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

  sf::RenderWindow &window;
  std::unordered_set<u64> set;
  std::vector<std::vector<double>> dist;

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
  auto &window = init_window(1.5);

  vrp_data_storage vrp(window);

  vrp.read_data("..\\input\\C108.txt");
  vrp.generate_initial_solution(vrp_data_storage::Heuristics::Greedy);

  auto regular_map = get_map("..\\input");
  auto bonus_map = get_map("..\\bonus");
  std::reference_wrapper<path_map_t> current_map = regular_map;

  bool bonus = false;
  sf::Event event;
  while (window.isOpen()) {
    vrp.draw();

    while (window.pollEvent(event)) {
      switch (event.type) {

      case sf::Event::KeyPressed:
        if (event.key.code == sf::Keyboard::Escape) { //  esc = close
          window.close();
        }
        if (event.key.code == sf::Keyboard::Space) { // change colors
          vrp.draw(true);
        } else if (event.key.code == sf::Keyboard::O) { // optimize
          vrp.perform_local_search();
          std::cout << "done" << std::endl;
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

          std::string name = std::filesystem::path(path).stem().string() +
                             ": " + std::to_string(vrp.overall_distance());

          window.setTitle(name);
          std::cout << "#" << n << ": " << name << std::endl;
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