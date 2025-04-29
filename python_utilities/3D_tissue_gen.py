

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from random import randrange, uniform
import skimage
import os
from skimage import data, color, transform
from skimage import exposure
from skimage.measure import label, regionprops, regionprops_table
from skimage import filters, morphology, measure


class Cell:
    def __init__(self, radius, position,cell_array):
        self.radius = radius
        self.position = np.array([position[0],position[1]])
        self.overlap = False
        self.area=0
        self.cell_array=cell_array
        self.volume=0
    def get_area(self):
        self.area=3.14159*self.radius ** 2
        return self.area
    def get_volume(self):
        self.volume=(4/3)*3.14159*self.radius ** 3
        return self.volume
    def is_overlap(self):
        for cell in self.cell_array:
            extra_spacing=0
            spacing=max(self.radius,cell.radius)*uniform(0,0.5)-extra_spacing
            if np.linalg.norm(self.position-cell.position) < (self.radius+cell.radius-spacing):
                self.overlap =True
            elif np.linalg.norm(self.position-cell.position)< self.radius+cell.radius and spacing>0:
                self.area=self.area-intersection_of_two_circles(self.position,self.radius,cell.position,cell.radius)
        return self.overlap
    def add_to_list(self):
        # circle_elements=[np.array([self.position[0],self.position[1]]),self.radius,self.area]
        self.cell_array.append(self)
        return self.cell_array

class Fiber:
    def __init__(self, radius, length, position, fiber_array,angle=0):
        self.radius = radius
        self.length = length
        self.position = position
        self.overlap = False
        self.area=0
        self.fiber_array=fiber_array

        self.x = position[0]
        self.y = position[1]
        self.width = 2*radius
        self.height = length
        self.angle = angle
        self.volume=0
    def get_area(self):
        self.area=self.width*self.height
        return self.area
    def get_volume(self):
        self.volume=self.length*3.1159*self.radius**2
        return self.volume
    # def is_overlap(self):
    #     for i in range(len(self.cell_array)):
    #         spacing=max(self.radius,self.cell_array[i][1])*uniform(0,0.5)
    #         if np.linalg.norm(self.position-self.cell_array[i][0]) < (self.radius+self.cell_array[i][1]-spacing):
    #             self.overlap =True
    #         elif np.linalg.norm(self.position-self.cell_array[i][0])< self.radius+self.cell_array[i][1] and spacing>0:
    #             self.area=self.area-intersection_of_two_circles(self.position,self.radius,self.cell_array[i][0],self.cell_array[i][1])
    #     return self.overlap
    #

    def get_corners(self):
        """
        Get the coordinates of the rectangle's corners after rotation.
        
        :return: List of (x, y) tuples representing the corners
        """
        half_width = self.width / 2
        half_height = self.height / 2
        
        # Define the unrotated corners relative to the center
        corners = [
            (-half_width, -half_height),
            (half_width, -half_height),
            (half_width, half_height),
            (-half_width, half_height)
        ]
        
        # Rotate the corners
        angle_rad = math.radians(self.angle)
        rotated_corners = []
        for (dx, dy) in corners:
            x_rot = dx * math.cos(angle_rad) - dy * math.sin(angle_rad)
            y_rot = dx * math.sin(angle_rad) + dy * math.cos(angle_rad)
            rotated_corners.append((self.x + x_rot, self.y + y_rot))
        
        return rotated_corners

    def overlaps_rectangle(self, other):
        """
        Check if this rectangle overlaps with another rectangle.
        
        :param other: Another Cylinder2D object
        :return: True if overlapping, False otherwise
        """
        # Use Separating Axis Theorem (SAT) for rectangle-rectangle collision detection
        def project(corners, axis):
            dots = [x * axis[0] + y * axis[1] for (x, y) in corners]
            return min(dots), max(dots)

        def overlap(proj1, proj2):
            return not (proj1[1] < proj2[0] or proj2[1] < proj1[0])

        # Get the corners of both rectangles
        corners1 = self.get_corners()
        corners2 = other.get_corners()

        # Define the axes to test (normals of the edges)
        axes = []
        for i in range(len(corners1)):
            x1, y1 = corners1[i]
            x2, y2 = corners1[(i + 1) % len(corners1)]
            edge = (x2 - x1, y2 - y1)
            normal = (-edge[1], edge[0])
            length = math.hypot(normal[0], normal[1])
            if length > 0:
                axes.append((normal[0] / length, normal[1] / length))

        for i in range(len(corners2)):
            x1, y1 = corners2[i]
            x2, y2 = corners2[(i + 1) % len(corners2)]
            edge = (x2 - x1, y2 - y1)
            normal = (-edge[1], edge[0])
            length = math.hypot(normal[0], normal[1])
            if length > 0:
                axes.append((normal[0] / length, normal[1] / length))

        # Check for overlap on all axes
        for axis in axes:
            proj1 = project(corners1, axis)
            proj2 = project(corners2, axis)
            if not overlap(proj1, proj2):
                return False
        return True

    def overlaps_cell(self, cell):
        """
        Check if this rectangle overlaps with a circle.
        
        :param circle: A tuple (x, y, radius) representing the circle
        :return: True if overlapping, False otherwise
        """
        circle_x, circle_y, radius = cell
        
        # Transform circle coordinates into rectangle's local space
        angle_rad = -math.radians(self.angle)
        dx = circle_x - self.x
        dy = circle_y - self.y
        local_x = dx * math.cos(angle_rad) - dy * math.sin(angle_rad)
        local_y = dx * math.sin(angle_rad) + dy * math.cos(angle_rad)
        
        # Find the closest point on the rectangle to the circle
        closest_x = max(-self.width / 2, min(local_x, self.width / 2))
        closest_y = max(-self.height / 2, min(local_y, self.height / 2))
        
        # Calculate the distance between the circle and the closest point
        distance_sq = (local_x - closest_x) ** 2 + (local_y - closest_y) ** 2
        return distance_sq <= radius ** 2

    # def draw(self, ax, color="blue"):
    #     """
    #     Draw the rectangle using matplotlib.patches.
    #
    #     :param ax: Matplotlib axis object
    #     :param color: Color of the rectangle (default: blue)
    #     """
    #     rect = patches.Rectangle(
    #         (self.x - self.width / 2, self.y - self.height / 2),
    #         self.width, self.height, angle=self.angle,
    #         color=color, fill=True, alpha=0.5
    #     )
    #     ax.add_patch(rect)
    def add_to_list(self):
        # fiber_elements=[np.array([self.position[0],self.position[1]]),self.radius, self.area, self.length, self.angle ]
        # self.fiber_array.insert(-1,fiber_elements)
        self.fiber_array.append(self)
        return self.fiber_array


def intersection_of_two_circles(position1,radius1,position2,radius2):
    d=np.linalg.norm(position2-position1)
    r=radius2
    R=radius1
    Area=(r**2)*math.acos(((d**2)+(r**2)-(R**2))/(2*d*r))+(R**2)*math.acos(((d**2)+(R**2)-(r**2))/(2*d*R))-0.5*((-1*d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))**(1/2)
    return Area 

def make_random_cells(cell_array,total_area,radius_upper_bound,radius_lower_bound):
    x_lower_bound=-300;
    x_upper_bound=300;

    y_lower_bound=-150;
    y_upper_bound=150;

    x_rand = uniform(x_lower_bound,x_upper_bound)
    y_rand = uniform(y_lower_bound,y_upper_bound)
    radius_rand= uniform(radius_lower_bound,radius_upper_bound)
    position=(x_rand,y_rand)
    a=Cell(radius_rand,position,cell_array) ##make a random circle  
    a.get_area() ##calculate the circle area
    a.is_overlap() ##check if it overlaps another circle by the amount defined in the circle class
    if a.overlap==False and (x_rand-radius_rand)>=x_lower_bound and (x_rand+radius_rand)<=x_upper_bound:## add check that it is not overlaping or on the boundary
        if (y_rand-radius_rand)>=y_lower_bound and (y_rand+radius_rand)<=y_upper_bound:
            a.add_to_list()
            total_area=total_area+a.area
    return total_area

def make_random_cells_and_fibers(cell_array,fiber_array,total_area,radius_upper_bound,radius_lower_bound,fiber_radius_lower_bound,fiber_radius_upper_bound):
    x_lower_bound=-300;
    x_upper_bound=300;

    y_lower_bound=-150;
    y_upper_bound=150
    # x_lower_bound=-150
    # x_upper_bound=150
    #
    # y_lower_bound=-75
    # y_upper_bound=75
    length =12
    for i in range(1):
        x_rand = uniform(x_lower_bound,x_upper_bound)
        y_rand = uniform(y_lower_bound,y_upper_bound)
        radius_rand= uniform(fiber_radius_lower_bound,fiber_radius_upper_bound)
        position=(x_rand,y_rand)
        rand_angle=uniform(0,360)
    #def __init__(self, radius, length, position, angle=0, cell_array):
        b=Fiber(radius_rand,length,position,fiber_array,rand_angle)#make a random fiber
        no_cell_overlap=True
        no_fiber_overlap=True
        inside_bound=True
        for cell in cell_array:
            cell_param=(cell.position[0],cell.position[1],cell.radius) #position and radius
            if b.overlaps_cell(cell_param):
                no_cell_overlap=False
        for fiber in fiber_array:
            if b.overlaps_rectangle(fiber):
                no_fiber_overlap=False
        if no_cell_overlap!=False and no_fiber_overlap!=False:
            for corner in b.get_corners():
                 corner=corner+position
                 if corner[0] < x_lower_bound or corner[0] > x_upper_bound:
                    inside_bound=False
                 if corner[1] < y_lower_bound or corner[1] > y_upper_bound:
                    inside_bound=False
            if inside_bound:
                b.get_area()
                b.add_to_list()
                total_area=total_area+b.area
    x_rand = uniform(x_lower_bound,x_upper_bound)
    y_rand = uniform(y_lower_bound,y_upper_bound)
    radius_rand= uniform(radius_lower_bound,radius_upper_bound)
    position=(x_rand,y_rand)
    a=Cell(radius_rand,position,cell_array) ##make a random circle
    cell_param=(a.position[0],a.position[1],a.radius) #position and radius
    a.get_area() ##calculate the circle area
    a.is_overlap() ##check if it overlaps another circle by the amount defined in the circle class
    fiber_overlap =False
    for fiber in fiber_array:
        if fiber.overlaps_cell(cell_param):
            no_fiber_overlap=False
    if a.overlap==False and no_fiber_overlap!=False and (x_rand-radius_rand)>=x_lower_bound and (x_rand+radius_rand)<=x_upper_bound:## add check that it is not overlaping or on the boundary
        if (y_rand-radius_rand)>=y_lower_bound and (y_rand+radius_rand)<=y_upper_bound:
            a.add_to_list()
            total_area=total_area+a.area
    return total_area


def plot2(ax, list_of_circles,artists):
    for cell in list_of_circles:
        circ = mpatches.Circle(cell.position, color="green", ec="black")
        ax.add_patch(circ)

def plot_fibers(ax, list_of_fibers,artists):
    for fiber in list_of_fibers:
        rect = mpatches.Rectangle(
            (fiber.x - fiber.width / 2, fiber.y - fiber.height / 2),
            fiber.width, fiber.height, angle=fiber.angle,
            color="black", fill=True, alpha=1
        )
        ax.add_patch(rect)
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import math
#
# class Cylinder2D:
#     def __init__(self, x, y, width, height, angle=0):
#         """
#         Initialize a 2D rectangle (representing a cylinder).
#
#         :param x: x-coordinate of the rectangle's center
#         :param y: y-coordinate of the rectangle's center
#         :param width: width of the rectangle
#         :param height: height of the rectangle
#         :param angle: rotation angle of the rectangle in degrees (default: 0)
#         """
#         self.x = x
#         self.y = y
#         self.width = width
#         self.height = height
#         self.angle = angle
#
#
# # Example usage
# if __name__ == "__main__":
#     # Create a figure and axis
#     fig, ax = plt.subplots()
#
#     # Create two Cylinder2D objects
#     rect1 = Cylinder2D(2, 2, 4, 2, angle=30)
#     rect2 = Cylinder2D(3, 3, 3, 3, angle=45)
#
#     # Create a circle (x, y, radius)
#     circle = (4, 4, 1.5)
#
#     # Draw the objects
#     rect1.draw(ax, color="blue")
#     rect2.draw(ax, color="green")
#     circle_patch = patches.Circle((circle[0], circle[1]), circle[2], color="red", alpha=0.5)
#     ax.add_patch(circle_patch)
#
#     # Check for overlaps
#     print("Rectangle 1 overlaps Rectangle 2:", rect1.overlaps_rectangle(rect2))
#     print("Rectangle 1 overlaps Circle:", rect1.overlaps_circle(circle))
#     print("Rectangle 2 overlaps Circle:", rect2.overlaps_circle(circle))
#
#     # Set axis limits and display
#     ax.set_xlim(0, 6)
#     ax.set_ylim(0, 6)
#     ax.set_aspect("equal")
#     plt.show()
def generate_csv(cell_list, fiber_list):
    with open("cell_list.csv", "w") as f:
        label_string=f'x,y,z,type,volume,custom:length,custom:angle'
        label_string=label_string+"\n"
        f.write(label_string)
        for cell in cell_list:
            volume=cell.get_volume()
            cell_string=f"{cell.position[0]},{cell.position[1]},0,granulosa,{volume},1,0"
            cell_string=cell_string+"\n"
            f.write(cell_string)
        for fiber in fiber_list:
            volume=fiber.get_volume()
            fiber_string=f'{fiber.position[0]},{fiber.position[1]},0,matrix,{volume},{fiber.length},{fiber.angle}'
            fiber_string=fiber_string+"\n"
            f.write(fiber_string)
def main():
    # test_nodes=[(0,-0.5),(0,-1),(0,-1.5),(-2,1),(0,-2),(0,-2.5),(0,-3),(0,-3.5),(0,-4),(0,0),(0,0.5),(0,1),(0,1.5),(-2,1),(0,2),(0,2.5),(0,3),(0,3.5),(0,4)]
    # test_centroid=[(-1,1.5)]
    # test_circle_center=[(8,1.5)]
    # test_radius=2
    # is_interior(test_nodes,test_centroid,test_circle_center,test_radius)
    list_of_circles=[]
    list_of_fibers=[]
    artists=[]
    total_area=0
    count=0
    x_lower_bound=-300;
    x_upper_bound=300;

    y_lower_bound=-150;
    y_upper_bound=150;
    side1=math.fabs(x_upper_bound)+math.fabs(x_lower_bound)
    side2=math.fabs(y_upper_bound)+math.fabs(y_lower_bound)
    remaining_area=side1*side2-total_area
    packing_density=total_area/(side1*side2)
    radius_lower_bound=4
    radius_upper_bound=7
    fiber_radius_lower_bound=1
    fiber_radius_upper_bound=2

    cell_max_rad=7
    test_points=[]
    while packing_density<0.6 and len(list_of_circles)<1600 and len(list_of_fibers)<4000 and count<9000:
        # total_area=make_random_circle(list_of_circles,total_area,1,4,c_hole)
        # total_area=make_random_cells(list_of_circles,total_area,4,7)
        total_area=make_random_cells_and_fibers(list_of_circles,list_of_fibers,total_area,radius_upper_bound,radius_lower_bound,fiber_radius_lower_bound,fiber_radius_upper_bound)
        count=count+1
        remaining_area=side1*side2-total_area
        packing_density=total_area/(side1*side2)
    print("free area",remaining_area)
    print("packing density", packing_density)
    # print(list_of_circles)
    generate_csv(list_of_circles,list_of_fibers)
    fig, ax = plt.subplots()
    plt.xlim(-500,500)
    plt.ylim(-500,500)
    ax.set_aspect('equal')
    title="packing cells with density of "+str(packing_density)
    ax.set_title(title)
    plot2(ax,list_of_circles,artists)
    plot_fibers(ax,list_of_fibers,artists)
    # plot(ax,list_of_circles,test_points,artists)
    for i in range(len(artists)):
        print(artists[i])
        ax.add_patch(artists[i])
    plt.show()
    print("\n")
    # print(frand)
# randrange gives you an integral value
#irand = randrange(0, 10)

# uniform gives you a floating-point value
#frand = uniform(0, 10)

if __name__ == "__main__":
    main()


