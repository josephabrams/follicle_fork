import numpy as np
import matplotlib.pyplot as plt 

start=0
stop=10
steps=20
x= np.linspace(start,stop,steps)
y= np.linspace(start,stop,steps)

y_test_1= np.sin(x)
y_test_2= x**2
###
dt=(stop-start)/steps
print(dt)
def calculate_new_position(x_current, x_old, acceleration, dt):
    return 2*x_current-x_old+acceleration*dt**2
def calculate_acceleration_test_1(position):
    return -1*np.sin(position)
def calculate_acceleration_test_2(position):
    return 2
def main():
    i=0
    x_old=y_test_1[0]
    x_current=y_test_1[1]
    for pos in x:
        # print(pos)
        acceleration=calculate_acceleration_test_1(pos)
        new_postion=calculate_new_position(x_current,x_old,acceleration,dt)
        # print(new_postion) 
        y[i]=new_postion
        x_old=x_current
        x_current=y[i]
        i=i+1

    fig, ax = plt.subplots()

    ax.plot(x, y+2.5, linewidth=2.0) 
    ax.plot(x, y_test_1, 'o-', linewidth=2)
    ax.set(xlim=(0, 10), xticks=np.arange(0, 10),
       ylim=(-5, 5), yticks=np.arange(-5, 5))
    plt.show()

if __name__ == "__main__":
    main()
