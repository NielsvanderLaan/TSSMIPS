#include "run.h"

bool find(int argc, char *argv[], string cmp)
{
  for (size_t idx = 1; idx != argc; ++idx)
  {
    string arg(argv[idx]);
    if (arg == cmp)
      return true;
  }
  return false;
}

bool use_rcuts(int argc, char *argv[])
{
  return not find(argc, argv, "OFF");
}

bool use_fenchel(int argc, char *argv[])
{
  return not find(argc, argv, "GOMORY");
}

bool solve_DEF(int argc, char *argv[])
{
  return find(argc, argv, "DEF");
}

bool solve_root(int argc, char *argv[])
{
  return find(argc, argv, "ROOT");
}

bool solve_tree(int argc, char *argv[])
{
  return find(argc, argv, "TREE");
}


int get_max_rounds(int argc, char *argv[])
{
  for (size_t idx = 1; idx != argc; ++idx)
  {
    string arg(argv[idx]);
    string mr("MAXROUNDS=");
    if (arg.find(mr) == 0)
      return stoi(arg.substr(mr.size(), arg.size() - mr.size()));
  }
  return -1;
}

int nthreads(int argc, char *argv[])
{
  int count = 1;
  for (size_t idx = 1; idx != argc; ++idx)
  {
    string arg(argv[idx]);
    string mr("NTHREADS=");
    if (arg.find(mr) == 0)
      count = stoi(arg.substr(mr.size(), arg.size() - mr.size()));
  }
  omp_set_num_threads(count);
  return count;
}



double get_time_limit(int argc, char *argv[])
{
  return 3 * 3600;
}

void details(vector<Type> types, int max_rounds, bool rcuts, bool fenchel, double time_limit, int thread_count)
{
  for_each(types.begin(), types.end(), [](Type type){cout << name(type) << "s\n";});
  cout << "reverse cuts: " << (rcuts ? "yes\n" : "no\n");
  cout << (fenchel ? "Fenchel" : "Gomory") << " cuts\n";
  cout << "Max rounds: " << max_rounds << '\n';
  cout << "Time limit: " << time_limit << "s\n";
  cout << "number of threads: " << thread_count << '\n';
  cout << endl;
}