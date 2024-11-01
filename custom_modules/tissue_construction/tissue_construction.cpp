#include "./tissue_construction.h"
#include <cmath>
//#define _USE_MATH_DEFINES

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
std::vector<std::vector<double>> create_spherical_shell(double cell_radius, double sphere_radius,double inner_radius) 
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


std::vector<std::vector<double>> twoD_symmetric_test_cells(double ring_radius) 
{
  std::vector<std::vector<double>> cells;
  std::vector<double> tempPoint(3, 0.0);
  double full_circle =2*M_PI;
  double angle_gap= M_PI/4;
  for (double angle = 0; angle < full_circle; angle += angle_gap) {
    tempPoint[0] = ring_radius*std::cos(angle);
    tempPoint[1] = ring_radius*std::sin(angle);
            cells.push_back(tempPoint);
  }
  return cells;
}
std::vector<std::vector<double>> x_test_rod(double cell_radius,std::vector<double> start_point,double length) 
{
  std::vector<std::vector<double>> cells;
  std::vector<double> tempPoint(3, 0.0);

  double x_spacing = 2*cell_radius;
  for (double x = start_point[0]; x < (start_point[0]+length); x += x_spacing) {
    tempPoint[0] = x;
    tempPoint[1] = start_point[1];
    tempPoint[2] = start_point[2];
            cells.push_back(tempPoint);
  }
  return cells;
}
//made the following test examples so that it's easy to modify them for quick tests
std::vector<std::vector<double>> two_cells(double cell_radius) 
{
  std::vector<std::vector<double>> cells;
  std::vector<double> tempPoint(3, 0.0);
  double x_spacing = 2*cell_radius;
  std::vector<double> point_1{-cell_radius,cell_radius,0.0};
  std::vector<double> x_length{x_spacing, 0.0, 0.0};
  std::vector<double> point_2=point_1+x_length;
  cells.push_back(point_1);
  cells.push_back(point_2);
  return cells;
}
std::vector<std::vector<double>> four_cells(double cell_radius) 
{
  std::vector<std::vector<double>> cells;
  std::vector<double> tempPoint(3, 0.0);
  double x_spacing = 2*cell_radius;
  double y_spacing = 2*cell_radius;
  std::vector<double> point_1{-cell_radius,-cell_radius,0.0};
  std::vector<double> x_length{x_spacing, 0.0, 0.0};
  std::vector<double> y_length{0.0, y_spacing, 0.0};
  std::vector<double> point_2=point_1+x_length;
  std::vector<double> point_3=point_1+y_length;
  std::vector<double> point_4=point_1+y_length+x_length;
  cells.push_back(point_1);
  cells.push_back(point_2);
  cells.push_back(point_3);
  cells.push_back(point_4);
  return cells;
}

std::vector<std::vector<double>> seven_cells(double cell_radius) 
{
  std::vector<std::vector<double>> cells;
  std::vector<double> tempPoint(3, 0.0);
  double x_spacing = 2*cell_radius;
  double y_spacing = 2*cell_radius;
  double z_spacing = 2*cell_radius;

  std::vector<double> center_cell{0.0,0.0,0.0};
  std::vector<double> x_length{x_spacing, 0.0, 0.0};
  std::vector<double> y_length{0.0, y_spacing, 0.0};
  std::vector<double> z_length{0.0, 0.0, z_spacing};
  std::vector<double> left=center_cell-x_length;
  std::vector<double> right=center_cell+x_length;
  std::vector<double> top=center_cell+y_length;
  std::vector<double> bottom=center_cell-y_length;
  std::vector<double> forward=center_cell+z_length;
  std::vector<double> back=center_cell-z_length;

  cells.push_back(center_cell);
  cells.push_back(left);
  cells.push_back(right);
  cells.push_back(top);
  cells.push_back(bottom);
  cells.push_back(forward);
  cells.push_back(back);
  return cells;
}
//Functions to check that cells are not passing into the oocyte or through the BM these are constraints on the allowed force
void outter_constraint(){

}
void inner_constraint(){

}
