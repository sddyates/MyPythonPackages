# coding: utf-8
import yt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import pylab
import matplotlib.pyplot as plt
from yt.utilities.physical_constants import mp, kb

def single_plot(ds):

    slc = yt.SlicePlot(ds, 'z' , field, center=center, width=width)
    #slc = yt.OffAxisSlicePlot(ds, L, field, center=center, north_vector=north_vector, width=width)
    """
    slc.annotate_streamlines('magnetic_field_x', 
                             'magnetic_field_z', 
                             density=2, 
                             plot_args={'color': 'black', 'linewidth': 1.0})
    """
    slc.set_cmap(field=field, cmap='jet')
    slc.set_zlim(field, 7.0e-17, 1.0e-20)
    slc.hide_colorbar()
    slc.hide_axes()
    #slc.annotate_grids(min_level=1, max_level=5)
    #slc.annotate_cell_edges(line_width=0.0001, alpha=0.5)

    slc.save('plots/highB_lowUV_closeup_topdown.png')

    
def time_series(ds):

    for ds in ts:
        slc = yt.SlicePlot(ds, 'z' , field, center=center, width=width)

        #slc = yt.ProjectionPlot(ds, 'z', field)
        #slc = yt.OffAxisSlicePlot(ds, L, field, center=center, north_vector=north_vector, width=width)

        slc.set_cmap(field=field, cmap='jet')
        """
        slc.annotate_streamlines('magnetic_field_x', 
                                 'magnetic_field_z', 
                                 density=2, 
                                 plot_args={'color': 'white', 'linewidth': 1.0})
        """
        """
        slc.annotate_streamlines('velocity_x', 
                                 'velocity_y', 
                                 density=2, 
                                 plot_args={'color': 'white', 'linewidth': 1.0})
        """                          
        #slc.annotate_line_integral_convolution('velocity_x', 'velocity_z', lim=(0.5,0.65))
        #slc.annotate_cell_edges()
        #slc.annotate_velocity()
        #slc.annotate_line((10.1, 0.0), (-15, 7.2), coord_system='plot')

        slc.save('plots/Free_AMR_higB_highUV/close_up_z.png')


def quad_plot(ds):

    fig = plt.figure()
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2, 2),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad="1.2%")

    for o in range(2):
        for i, fn in enumerate(fns):

            if i == 0 or i == 2:
                slc = yt.SlicePlot(ds, 
                                   'z' , 
                                   field, 
                                   center=center[o], 
                                   width=width[o])
                slc.annotate_streamlines('velocity_x', 
                                         'velocity_y', 
                                         density=1.0,
                                         factor=16, 
                                         plot_args={
                                             'color': 'white', 
                                             'linewidth': 0.25})
                slc.set_cmap(field=field, cmap='jet')
            else:
                slc = yt.OffAxisSlicePlot(ds, 
                                          L[o],
                                          field, 
                                          center=center[o], 
                                          north_vector=north_vector, 
                                          width=width[o])
                slc.annotate_streamlines('magnetic_field_x', 
                                         'magnetic_field_z', 
                                         density=1.0,
                                         factor=16, 
                                         plot_args={
                                             'color': 'white', 
                                             'linewidth': 0.25})
                slc.set_cmap(field=field, cmap='jet')

            if i != 2:
                slc.hide_axes()


            if o == 0:
                slc.set_zlim(field, 1e-20, 5e-15)
            else:
                slc.set_zlim(field, 1e-20, 7e-16)

            plot = slc.plots[field]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]

            slc._setup_plots()
        
        if o == 0:
            slc.save("plots/global.pdf")
        else:
            slc.save("plots/planet.pdf")


def double_plot(settings, ds):

    fig = plt.figure()
    grid = AxesGrid(fig, (0.09,0.09,0.8,0.8),
                    nrows_ncols = (1, 2),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad="1.2%")

    for i, fn in enumerate(fns):

        if settings["option"] == "top_down":
            slc = yt.SlicePlot(ds, 
                               'z' , 
                               settings["field"], 
                               center=settings["center"], 
                               width=settings["width"],
                               fontsize=settings["font"])
            if settings["streamlines"]:
                slc.annotate_streamlines('velocity_x', 
                                         'velocity_y', 
                                         density=1.5,
                                         factor=16, 
                                         plot_args={
                                             'color': 'black', 
                                             'linewidth': 0.25})

        if settings["option"] == "side_on":
            slc = yt.OffAxisSlicePlot(ds, 
                                      settings["L"],
                                      settings["field"], 
                                      center=settings["center"], 
                                      north_vector=settings["north_vector"], 
                                      width=settings["width"],
                                      fontsize=settings["font"])
            slc.set_xlabel('x $\ (\mathrm{R}_{\odot})$')
            slc.set_ylabel('z $\ (\mathrm{R}_{\odot})$')
            if settings["streamlines"]:
                slc.annotate_streamlines('magnetic_field_x', 
                                         'magnetic_field_z', 
                                         density=1.5, 
                                         factor=16, 
                                         plot_args={
                                             'color': 'white', 
                                             'linewidth': 0.75})
            
        slc.set_cmap(field=settings["field"], cmap='jet')
        slc.set_zlim(settings["field"], settings["lim"][0], settings["lim"][1])

        plot = slc.plots[settings["field"]]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        slc._setup_plots()

    slc.save("plots/"+settings["save_name"]+".pdf")



units_override = {"length_unit":(1,"Rsun"),
                  "time_unit":(6.955e+05 ,"s"),
                  "mass_unit":(3.36427433875e+17,"g"),
                  "magnetic_unit":(1.121e-02,"G")}

settings = {"option":"side_on",
            "field":"density",
            'north_vector':([0, 0, 1]),
            "center":[0, 0, 0],
            "width":((31.9, "Rsun"), (31.9, "Rsun")),
            "L":[0, -1, 0],
            "lim":[1.0e-21, 5.0e-15],
            "streamlines":True,
            "save_name":"side_on_planet",
            "font":14}


