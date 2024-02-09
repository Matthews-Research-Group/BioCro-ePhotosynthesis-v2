#include <iomanip>
#include <iostream>
#include <vector>

#include <nlopt.hpp>
int count=0;
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    ++count;
    std::cout << "this is the evalution # "<<count<<std::endl;  
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

int main()
{
  nlopt::opt opt(nlopt::AUGLAG, 2);
  nlopt::opt myopt(nlopt::LD_MMA, 2);
  std::vector<double> lb{-HUGE_VAL,1e-3};

  myopt.set_xtol_rel(1e-4);
  opt.set_local_optimizer(myopt);

  opt.set_lower_bounds(lb);
  opt.set_min_objective(myfunc, NULL);
  my_constraint_data data[2] = { {2,0}, {-1,1} };
  opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
  opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
  opt.set_xtol_rel(1e-4);
  opt.set_maxeval(10000);
  std::vector<double> x(2);
  x[0] = 1.234; x[1] = 5.678;
  double minf;

  try{
      nlopt::result result = opt.optimize(x, minf);
      std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
          << std::setprecision(10) << minf << std::endl;
  }
  catch(std::exception &e) {
      std::cout << "nlopt failed: " << e.what() << std::endl;
  }
  return 0;
}
