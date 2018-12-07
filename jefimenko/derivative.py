import numpy as np
import pdb

def currents_time_diff(grid):
    # this will calculate the time derivative of the currents
    # first make a list of all of the currents
    for n in range(len(grid.currents[0])):
        current_time_list = []
        # notice that amps is included in this so it will not be needed later
        for time in range(grid.time_size):
            current_time_list.append(grid.currents[time][n].direction
                                     * grid.currents[time][n].amps)
        # now calculate the gradient with numpy
        gradient = np.gradient(current_time_list, grid.delta_t, axis=0)

        # now put all of the values into grid.currents for storage
        for time in range(grid.time_size):
            grid.currents[time][n].diff_t = gradient[time]
