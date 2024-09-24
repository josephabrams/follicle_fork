import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def hookes_law(spring_k,rest_length, position_1=np.ndarray(shape=(1,3)), position_2=np.ndarray(shape=(1,3))):
    displacement=position_2-position_1
    length=np.sqrt(displacement[0]**2+displacement[1]**2+displacement[2]**2)
    unit_vec=displacement/np.sqrt(np.vdot(displacement,displacement))
    
    return spring_k*(length-rest_length)*unit_vec
    
def plot_cells(cell_1_position, cell_2_position, plot_name=' ', cell_1_radius=1, cell_2_radius=1):
    circle_1_pos=(cell_1_position[0],cell_1_position[1])
    circle_2_pos=(cell_2_position[0],cell_2_position[1])
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
def main():
    position_A=np.array([5,0,0])  
    position_B=np.array([10,5,0])  
    spring_constant=2.0
    mass1=1.0
    mass2=1.0
    plot_cells(position_A,position_B, 'before')
    test_result=hookes_law(spring_constant, 10, position_A, position_B)
    test_result1=test_result/mass1
    test_result2=test_result/mass2
    plot_cells((position_A+test_result1),(position_B-test_result2), 'after')
    print("TEST RESULT: ", test_result)

if __name__ == "__main__":
    main()
