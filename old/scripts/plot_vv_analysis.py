import numpy as np
import matplotlib.pyplot as plt
import signac
import gsd.hoomd

def plot_vv_analysis(job,crosslink_type):
    #load the data
    data = np.genfromtxt(job.fn(f'contract_bonds_analysis/voronoi_volumes_{crosslink_type}.txt'))

    #create a histogram
    if crosslink_type == 'all':
        color = 'pink'
    elif crosslink_type == 'thiol':
        color = 'orange'
    elif crosslink_type == 'ene':
        color = 'blue'
    else:
        raise ValueError(f"Unknown crosslink type: {crosslink_type}")
    plt.hist(x=data, bins=100, density=True, alpha=0.5, color=color, label='Voronoi volumes')
    plt.xlabel('Voronoi volume')
    plt.ylabel('Density')

    #add a vertical line for the mean
    mean = np.mean(data)
    plt.axvline(mean, color='r', linestyle='dashed', linewidth=1)
    y_max = plt.ylim()[1]
    plt.text(mean*1.1, 0.95*y_max, 'Mean: {:.2f}'.format(mean), color='red')

    #Calculate and write the standard deviation
    std = np.std(data)
    plt.text(x=mean*1.1, y=0.85*y_max, s='Standard deviation: {:.2f}'.format(std), color='black')

    #add a vertical line for the expected volume
    frame = gsd.hoomd.open(job.fn('polymerize.gsd'))[-1]
    volume = frame.configuration.box[0]*frame.configuration.box[1]*frame.configuration.box[2]
    expected_volume = volume/len(data)
    plt.axvline(expected_volume, color='g', linestyle='dashed', linewidth=1)
    plt.text(expected_volume*1.1, 0.9*y_max, 'Expected: {:.2f}'.format(expected_volume), color='green')

    plt.title(f'Voronoi volume distribution for {crosslink_type} crosslinkers')

    #TODO: save the plots however you want. Might need to change some stuff, 
    #   either in here or in the iteration
    # plt.savefig(job.fn(f'contract_bonds_analysis/voronoi_volumes_{crosslink_type}.png'))
    plt.show() # only for testing
    plt.close()

    return mean, std, expected_volume

#import the data
project = signac.get_project()
jobs = project.find_jobs()

for job in jobs:
    if 'fbc867' not in job.id:
        continue

    # plot for all crosslinkers
    mean, std, expected = plot_vv_analysis(job,crosslink_type = 'all')

    # plot for only thiol crosslinkers
    mean, std, expected = plot_vv_analysis(job,crosslink_type = 'thiol')

    # plot for only ene crosslinkers
    mean, std, expected = plot_vv_analysis(job,crosslink_type = 'ene')

    plt.show()