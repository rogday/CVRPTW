#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <list>
#include <unordered_set>

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

public:
  vrp_data_storage(std::string const &filename) {
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

    // std::cout << customers;

    generate_distance_matrix();
    generate_initial_solution();
  }

  void generate_initial_solution() {
    std::unordered_set<u64> set;
    for (u64 i = 0; i < customers.size(); ++i)
      set.insert(i);

    double overall_distance = 0.0;

    for (auto &path : paths) {

      i64 current_time = 0, warehouse_closing = customers[0].due_date;
      i64 capacity = vehicle_capacity;
      path.push_back(0);

      while (current_time < warehouse_closing) {
        i64 next = -1;
        double max_profit = 0.0;

        for (u64 i : set) {
          i64 open_time = customers[i].ready_time - current_time;

          double profit =
              profitability(path.back(), i, open_time,
                            customers[i].due_date - current_time, capacity);

          // std::cout << profit << std::endl;
          if (profit > max_profit) {
            max_profit = profit;
            next = i;
          }
        }

        if (next != -1) {

          current_time += i64(distance(path.back(), next) + 0.5);

          if (next != 0) {
            if (capacity >= customers[next].demand)
              set.erase(next);

            capacity -= std::min(capacity, customers[next].demand);
          } else
            capacity = vehicle_capacity;

          std::cout << next << " ";
          /*std::cout << ", time: " << current_time
                    << ", distance: " << distance(path.back(), next)
                    << ", capacity: " << capacity << ", i=" << path.back()
                    << ", j=" << next << std::endl;*/

          path.push_back(next);
        } else
          break;
      }

      overall_distance += current_time;
      std::cout << std::endl << std::endl;
    }
    std::cout << overall_distance << std::endl;
  }

private:
  double distance(u64 i, u64 j) {
    return std::sqrt(std::pow(customers[i].x - customers[j].x, 2) +
                     std::pow(customers[i].y - customers[j].y, 2));
  }

  // ->> max
  double profitability(u64 start, u64 end, i64 open_time, i64 close_time,
                       i64 capacity) {
    double d = distance(start, end);
    double ret = 0.0;
    open_time = (open_time <= 0) ? 1 : open_time;

    if (d == 0)
      return 0;

    if (close_time - d > 0) {
      ret = 1.0 / (open_time) + max_distance / d;

      if (capacity > customers[end].demand)
        ret += 1.0;
      else {
        if (customers[end].demand != 0)
          ret += capacity / customers[end].demand;
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
  std::vector<std::vector<double>> distance_matrix;

  std::vector<std::list<u64>> paths;
};

int main() {
  vrp_data_storage vrp("..\\input\\C108.txt");

  return 0;
}
