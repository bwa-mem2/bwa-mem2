#include <vector>
#include <iostream>
#include <fstream>
#include "rmi.h"

int main() {
  // load the data
  std::vector<uint64_t> data;
  std::ifstream in("../osm_cellids_200M_uint64",
                   std::ios::binary);
  
  // Read size.
  uint64_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
  data.resize(size);
  // Read values.
  in.read(reinterpret_cast<char*>(data.data()), size*sizeof(uint64_t));
  in.close();

  std::cout << "Data loaded." << std::endl;

  std::cout << "RMI status: " << rmi::load("rmi_data") << std::endl;

  size_t err;
  
  for (uint64_t key_index = 0; key_index < size; key_index++) {
    uint64_t lookup = data[key_index];
    uint64_t true_index = (uint64_t)
      std::distance(data.begin(), std::lower_bound(data.begin(),
                                                   data.end(),
                                                   lookup));
    uint64_t rmi_guess = rmi::lookup(lookup, &err);
    
    uint64_t diff = (rmi_guess > true_index ? rmi_guess - true_index : true_index - rmi_guess);
    if (diff > err) {
      std::cout << "Search key: " << lookup
                << " Key at " << true_index << ": " << data[true_index] 
                << " RMI guess: " << rmi_guess << " +/- " << err
                << " diff: " << diff << std::endl;
      exit(-1);
    }
  }
  
  rmi::cleanup();
  exit(0);
}
