
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class Cell:
  def __init__(self, radius, position,old_position, mass, youngs, poisson):
    self.radius = radius
    self.position = position
    self.old_position= old_position
    self.mass = mass
    self.youngs=youngs
    self.poisson=poisson


def youngs_modulus(E1,E2,poisson1,poisson2,radius1,radius2, position_1=np.ndarray(shape=(1,3)), position_2=np.ndarray(shape=(1,3))):
    displacement=position_2-position_1
    penetration_depth= (displacement[0,0]**2+displacement[0,1]**2+displacement[0,2]**2)**(1/2)-radius1-radius2
    print(displacement[0,0]**2)
    print(penetration_depth)
    if penetration_depth>0 :
        return np.array([0,0,0])
    else:
        constant_term= ((3*math.pi)**(2/3))/2
        unit_vec=displacement/np.sqrt(np.vdot(displacement,displacement))
        V1=(1-poisson1**2)/(math.pi*E1)
        V2=(1-poisson2**2)/(math.pi*E2)
        D1= 1/radius1*2 
        D2= 1/radius2*2 
        Vterm= (V1+V2)**(2/3)
        Dterm= (D1+D2)**(1/3)
        force_magnitude=1*penetration_depth/(constant_term*Vterm*Dterm)
        unit_vec=displacement/np.sqrt(np.vdot(displacement,displacement))
        return force_magnitude*unit_vec
def cell_youngs_modulus(Cell1,Cell2):
    return youngs_modulus(Cell1.youngs,Cell2.youngs,Cell1.poisson, Cell2.poisson,Cell1.radius,  
                          Cell2.radius,Cell1.position,Cell2.position)

def plot_cells(cell_1_position, cell_2_position, plot_name=' ', cell_1_radius=1, cell_2_radius=1):
    circle_1_pos=(cell_1_position[0,0],cell_1_position[0,1])
    circle_2_pos=(cell_2_position[0,0],cell_2_position[0,1])
    xy_artists = [
        mpatches.Circle(circle_1_pos, radius=cell_1_radius,alpha=0.2, ec="black", fc='black'),
        mpatches.Circle(circle_2_pos, radius=cell_2_radius,alpha=1, ec="blue", fc='blue'),
    ]
    fig,ax=plt.subplots()
    for i in xy_artists:
        ax.add_patch(i)
    ax.autoscale_view()
    
    ax.set_aspect('equal', 'box')
    ax.set_title(plot_name)
    plt.show()
def plot_cells_and_save(cell_1_position, cell_2_position, plot_name=' ', cell_1_radius=1, cell_2_radius=1):
    circle_1_pos=(cell_1_position[0,0],cell_1_position[0,1])
    circle_2_pos=(cell_2_position[0,0],cell_2_position[0,1])
    xy_artists = [
        mpatches.Circle(circle_1_pos, radius=cell_1_radius,alpha=0.2, ec="black", fc='black'),
        mpatches.Circle(circle_2_pos, radius=cell_2_radius,alpha=1, ec="blue", fc='blue'),
    ]
    fig,ax=plt.subplots()
    for i in xy_artists:
        ax.add_patch(i)
    # ax.autoscale_view()
    ax.set_xlim([-20,20])
    ax.set_ylim([-20,20])
    ax.set_aspect('equal', 'box')
    ax.set_title(plot_name)
    plt.savefig(plot_name)
def calculate_new_position(x_current, x_old, acceleration, dt):
    return 2*x_current-x_old+acceleration*dt**2
def calc_first_position(x_current,acceleration,dt):
    return x_current+acceleration+(dt)**2 #forward euler twice
def acceleration_of_cells(dt,end_time, cell1,cell2):
    steps=int(end_time/dt)
    x=np.linspace(0,end_time,steps)
    step=0
    for time in x:
        force=cell_youngs_modulus(cell1,cell2)
        acceleration_1=force/cell1.mass
        acceleration_2=-1*force/cell2.mass
        print("FORCE: ", force)
        print("Acceleration: ", acceleration_1)
        # if time<dt :
        #     cell1.position=calc_first_position(cell1.old_position,acceleration_1,dt)
        #     cell2.position=calc_first_position(cell2.old_position,acceleration_2,dt)
        # else :
        cell1.position=calculate_new_position(cell1.position, cell1.old_position, acceleration_1, dt)
        # cell2.position=calculate_new_position(cell2.position, cell2.old_position, acceleration_2, dt)
        if step<10:
            plot_name="./output/plot00"+str(step)+".svg"
        elif step>=10 and step<100:
            plot_name="./output/plot0"+str(step)+".svg"
        else:
            plot_name="./output/plot"+str(step)+".svg"
        plot_cells_and_save(cell1.position, cell2.position, plot_name, cell1.radius, cell2.radius) 
        cell1.old_position=cell1.position
        # cell2.old_position=cell2.position
        step=step+1
        print("Cell1: ", cell1.old_position)
        print("Cell2: ", cell2.old_position)
        print("DONE MAKING IMAGES")
def main():

    position_A=np.array([[5,0,0]])  
    position_B=np.array([[10,5,0]])
    radius1=5
    radius2=5
    youngs1=385
    youngs2=385
    poisson1=0.5
    poisson2=0.5
    mass1=524
    mass2=524
    print(position_A,"shape: ", np.shape(position_A))
    print(position_B,"shape: ", np.shape(position_B))
    # plot_cells(position_A,position_B, 'before', cell_1_radius=radius1, cell_2_radius=radius2)
    test_result=youngs_modulus(youngs1,youngs2, poisson1,poisson2,radius1,radius2, position_A, position_B)
    test_result1=test_result/mass1
    test_result2=test_result/mass2
    # plot_cells((position_A+test_result1),(position_B-test_result2),'after: ',cell_1_radius=radius1,cell_2_radius=radius2)
    print("A (start) ", position_A)
    print("A (end) ", position_A+test_result1)
    print("B (start) ", position_B)
    print("B (end) ", position_B-test_result2)
  #def __init__(self, radius, position,old_position, mass, youngs, poisson):
    cell1= Cell(radius1,position_A,position_A, mass1,youngs1,poisson1)
    cell2= Cell(radius2,position_B,position_B, mass2,youngs2,poisson2)
    force=cell_youngs_modulus(cell1,cell2)
    acceleration_1=force/cell1.mass
    acceleration_2=-1*force/cell2.mass
    print("FORCE: ", force)
    print("Acceleration: ", acceleration_1)
    cell1.position=calc_first_position(cell1.old_position,acceleration_1,0.1)
    # cell2.position=calc_first_position(cell2.old_position,acceleration_2,0.1)
    acceleration_of_cells(0.1,100,cell1,cell2)

if __name__ == "__main__":
    main()
