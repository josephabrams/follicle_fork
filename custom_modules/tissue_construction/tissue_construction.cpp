#include "./tissue_construction.h"
std::vector<std::vector<double>> create_spheroid_2D(double cell_radius, double sphere_radius) 
{
  std::vector<std::vector<double>> cells;
  int xc = 0, yc = 0, zc = 0;
  double x_spacing = cell_radius * sqrt(3);
  double y_spacing = cell_radius * 2;
  std::vector<double> tempPoint(3, 0.0);

  
    for (double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++) {
      for (double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++) {
        tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
        tempPoint[1] = y + (xc % 2) * cell_radius;
        tempPoint[2] = 0;

        if (sqrt(norm_squared(tempPoint)) < sphere_radius) {
          /* output file of initial positions
          std::ofstream ofs;
          ofs.open ("temp_points_location.csv", std::ofstream::out |
          std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint)))<<"," <<"\n";
          ofs.close();
          */
          cells.push_back(tempPoint);
        }
      }
    }
  return cells;
}
std::vector<std::vector<double>> create_spheroid(double cell_radius, double sphere_radius) 
{
  std::vector<std::vector<double>> cells;
  int xc = 0, yc = 0, zc = 0;
  double x_spacing = cell_radius * sqrt(3);
  double y_spacing = cell_radius * 2;
  double z_spacing = cell_radius * sqrt(3);
  std::vector<double> tempPoint(3, 0.0);

  for (double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++) {
    for (double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++) {
      for (double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++) {
        tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
        tempPoint[1] = y + (xc % 2) * cell_radius;
        tempPoint[2] = z;

        if (sqrt(norm_squared(tempPoint)) < sphere_radius) {
          /* output file of initial positions
          std::ofstream ofs;
          ofs.open ("temp_points_location.csv", std::ofstream::out |
          std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint)))<<"," <<"\n";
          ofs.close();
          */
          cells.push_back(tempPoint);
        }
      }
    }
  }
  return cells;
}
std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius,double inner_radius) 
{
  std::vector<std::vector<double>> cells;
  int xc = 0, yc = 0, zc = 0;
  double x_spacing = cell_radius * sqrt(3);
  double y_spacing = cell_radius * 2;
  double z_spacing = cell_radius * sqrt(3);
  std::vector<double> tempPoint(3, 0.0);

  for (double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++) {
    for (double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++) {
      for (double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++) {
        tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
        tempPoint[1] = y + (xc % 2) * cell_radius;
        tempPoint[2] = z;

        if (sqrt(norm_squared(tempPoint)) < sphere_radius) {
          if (sqrt(norm_squared(tempPoint)) > inner_radius) {
            /* output file of initial positions
            std::ofstream ofs;
            ofs.open ("temp_points_location.csv", std::ofstream::out |
            std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint)))<<"," <<"\n";
            ofs.close();
            */
            cells.push_back(tempPoint);
          }
        }
      }
    }
  }
  return cells;
}
