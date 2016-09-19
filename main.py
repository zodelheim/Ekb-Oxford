
import numpy as np
from scipy.integrate import odeint
import os
from Ekb_Oxford_old_coop import Ekb_Oxford
from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.models.widgets import Slider, DataTable, TableColumn, PreText
from bokeh.io import curdoc, output_file
from bokeh.models import ColumnDataSource
import matplotlib.pyplot as plt


def plotting_1_param():

    var_names = ['v', 'd', 'f2', 'f2ds', 'f', 'w', 'ActFrac', 'ProdFrac', 'N', 'K_o', 'h', 'm', 'V_fibro',
                 'A', 'B_1', 'B_2', 'Ca_ds', 'Ca_i', 'Ca_rel', 'Ca_up', 'K_i', 'Na_i', 'l_1', 'l_2', 'l_3',
                 'V', 'xr1', 'xr2', 'xs', 'r', 's', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    directory = './Results/changes in s_c/'
    print 'hw'
    output_file("plots.html")
    A = np.load(directory+'data.npy')
    time = np.load(directory+'time.npy')
    changes = np.loadtxt(directory+'changes.txt')
    source = ColumnDataSource(data=dict(time=time, Y=A[0, 17, :]))
    Curve = Slider(start=0, end=len(A[0, :, 1])-1, value=17, step=1, title="plot", orientation = 'vertical')
    param = Slider(start=0, end=int(len(A[:,0,0])-1), step=1, value=0)
    p = figure(title= 'v - time', x_axis_label='time', y_axis_label='v', webgl=True)
    p.line('time', 'Y', source=source, line_width=3)
    data = dict(
        number = range(0, len(A[0, :, 1])),
        value = var_names
    )
    table_source = ColumnDataSource(data)
    table_columns = [
        TableColumn(field = 'number', title = 'number'),
        TableColumn(field = 'value', title = 'Curve')
    ]
    data_table = DataTable(source=table_source, columns=table_columns, width=400, height=600)
    # text
    text = PreText(text = 'A*'+str(changes[0]))
    def update_data(attr, old, new):
        Curve_number = int(Curve.value)
        param_number = int(param.value)
        source.data = dict(time=time, Y=A[param_number, Curve_number, :])
        p.title.text = var_names[Curve_number] + ' - time'
        p.yaxis.axis_label = str(var_names[Curve_number])
        param.title = Curve.title +'*'+str(changes[Curve_number])
        text.text =  'A*'+str(changes[param_number])
    Curve.on_change('value', update_data)
    param.on_change('value', update_data)
    l = layout([
        [p, Curve, data_table],
        [text, param]
    ])
    curdoc().add_root(l)
    curdoc().title = "Sliders"
    # show(l)

def plotting_2_param():
    directory = './'
    print 'hw'
    A = np.load(directory + 'data.npy')
    time = np.load(directory + 'time.npy')
    changes = np.loadtxt(directory + 'changes.txt')
    source = ColumnDataSource(data=dict(time=time, Y=A[0, 0, :]))
    Curve_on = Slider(start=0, end=len(A[0, :, 1]) - 1, value=0, step=1, title="param_on", orientation='vertical')
    Curve_off = Slider(start=0, end=int(len(A[:, 0, 0]) - 1), step=1, value=0, title="param _off")
    p = figure(title='Ca - time', x_axis_label='time', y_axis_label='v', webgl=True)
    p.line('time', 'Y', source=source, line_width=3)
    # text
    a_off = 200.0  # per_second (in intracellular_calcium_concentration)
    a_on= 70000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    b_1_off = 182.0  # per_second (in intracellular_calcium_concentration)
    b_1_on = 100000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    b_2_off = 3.0  # per_second (in intracellular_calcium_concentration)
    b_2_on = 1000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    on_value = a_on
    off_value = a_off
    text_on = PreText(text = 'value * '+str(changes[0]))
    text_off = PreText(text = 'value * '+str(changes[0]))
    # function to update data
    def update_data(attr, old, new):
        Curve_on_number = int(Curve_on.value)
        Curve_off_number = int(Curve_off.value)
        source.data = dict(time=time, Y=A[Curve_off_number, Curve_on_number, :])
        Curve_on.title = Curve_on.title + '*' + str(changes[Curve_on_number])
        Curve_off.title = Curve_off.title + '*' + str(changes[Curve_off_number])
        text_on.text = 'value * ' + str(changes[Curve_on_number])
        text_off.text = 'value * ' + str(changes[Curve_off_number])
    Curve_on.on_change('value', update_data)
    Curve_off.on_change('value', update_data)
    # inputs = widgetbox([Curve], [param])
    l = layout([
        [p, Curve_on],
        [Curve_off, text_off, text_on]
    ])
    curdoc().add_root(l)
    curdoc().title = "Sliders"
    # show(l)

def zero_dimentional():

    # # for old cooperativity
    Y_0 = np.array(
        [0.0, 0.0, 1.0, 0.997, 1.0, 0.0, 0.001914, 0.2854569, 7.455e-8, 3.988, 0.995, 0.0015, -20.0, 0.00015, 0.0, 0.0,
         2.55e-6, 6.15e-6, 0.989665, 0.994579, 139.054, 5.18787513189509, 0.436333342969918, 0.436333525334166,
         0.088805830771694, -93.658148, 8.88859784542779e-6, 1.53745791069154e-7, 0.001, 1.63117173173398e-8,
         0.997044616031121, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

    var_names = ['v', 'd', 'f2', 'f2ds', 'f', 'w', 'ActFrac', 'ProdFrac', 'N', 'K_o', 'h', 'm', 'V_fibro',
                 'A', 'B_1', 'B_2', 'Ca_ds', 'Ca_i', 'Ca_rel', 'Ca_up', 'K_i', 'Na_i', 'l_1', 'l_2', 'l_3',
                 'V', 'xr1', 'xr2', 'xs', 'r', 's']

    t_end = 10.
    t_start = 0.
    steps = 20000
    t = np.linspace(t_start, t_end, steps)
    A = np.zeros([46, steps])
    h = t[1] - t[0]
    a_off = 200.0  # per_second (in intracellular_calcium_concentration)
    a_on = 70000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    # a_off = 19.6  # per_second (in intracellular_calcium_concentration)
    # a_on = 32.7*1000  # per_millimolar_second (in intracellular_calcium_concentration)

    # b_1_off = 182.0  # per_second (in intracellular_calcium_concentration)
    # b_1_on = 100000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    # b_2_off = 3.0  # per_second (in intracellular_calcium_concentration)
    # b_2_on = 1000.0  # per_millimolar_second (in intracellular_calcium_concentration)
    # B_1_tot = 0.08  # millimolar (in intracellular_calcium_concentration)
    # B_2_tot = 0.1  # millimolar (in intracellular_calcium_concentration)

    # b_off = np.array([0.032, 3.33, 238., 0.46, 0.057, 60., 1300., 1300., 30., 30., 60., 60., 60., 110., 110., 110., 65000.])
    # b_on = np.array([2.37, 0.003, 34., 13.8, 0.0157, 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100.])*1000.
    # B_tot = np.array([140., 140., 24., 140., 140., 19., 37.4, 4.6, 13.4, 1.65, 25., 0.77, 0.02, 25., 0.77, 0.02, 140])/1000.

    # b_off = np.array([0.032, 3.33, 238., 0.46, 0.057, 60., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    # b_on = np.array([2.37, 0.003, 34., 13.8, 0.0157, 100., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])*1000.
    # B_tot = np.array([140., 140., 24., 140., 140., 19., 0.0, 0.0, 0.0, 0.0, 0., 0.0, 0.0, 0., 0.0, 0.0, 0.])/1000.
    b_off = np.array([182.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    b_on = np.array([100000.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    B_tot = np.array([0.08, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    k_A = 40.0  # per_millimolar (in parameters_izakov_et_al_1991)
    k_mu = 0.6  # dimensionless (in parameters)
    mu = 3.0    # dimensionless (in parameters)
    s_c = 35.0
    beta_3 = 0.015  # millinewton (in parameters)
    alpha_3 = 48.0  # per_micrometre (in parameters)
    var_array = np.logspace(-1, 1, 10)
    pars_number = 1
    if pars_number == 1:
        # Result_A = np.zeros([len(var_array), len(Y_0), steps])
        Result_A = np.zeros([len(var_array), steps/10])
        # for i in xrange(0, len(var_array)):
        # for i, var in enumerate(var_array):
        for i,var in enumerate([5]):
            print var
            # b_on[0] = 100000. * var_array[i]
            # b_off[0] = 182. * var_array[i]
            # B_tot[0] = 0.08 * var_array[i]
            # b_on[1] = 1000. * var_array[9]
            # b_off[1] = 3. * var_array[0]
            # B_tot[1] = 0.1 * var_array[i]
            A = odeint(Ekb_Oxford, Y_0, t,
                       args = (b_on, b_off,
                               B_tot,
                               var))
            Result_A[i, :] = [A[a, 17] for a in xrange(0, steps, 10)]

            F_XSE = beta_3 * (np.exp(alpha_3 * A[:, 24]) - 1.0)
            plt.figure(1)
            plt.plot(range(0, steps, 10), Result_A[i, :], label='par * ' + str(var_array[i]))
            plt.title('Ca_c - time')
            plt.legend()
            # plt.figure(2)
            # plt.plot(t, F_XSE, label='par * ' + str(var_array[i]))
            # plt.title('Force - time')
            # plt.legend()
    elif pars_number == 2:
        # # for 2 params change
        Result_A = np.zeros([len(var_array), len(var_array), steps/10])

        for i in xrange(len(var_array)):
            # for j in xrange(0, len(var_array)):
            for j in [1]:
                # b_on[0] = 100000. * var_array[j]
                # b_off[0] = 182. * var_array[j]
                # b_on[1] = 1000. * var_array[j]
                # b_off[1] = 3. * var_array[i]
                B_tot[0] = 0.08 * var_array[j]
                B_tot[1] = 0.1 * var_array[i]
                print j, i
                A = odeint(Ekb_Oxford, Y_0, t,
                           args = (b_on, b_off,
                                   B_tot,
                                   var_array[i]))

                # Result_A[i, j, :] = A[:, 17]
                Result_A[i, j, :] = [A[a, 17] for a in xrange(0, steps, 10)]
    ###############################################################################################
    np.savetxt('changes.txt', var_array)
    np.save('data', Result_A)
    np.save('time', t)
    plt.show()
    # F_XSE = np.zeros(steps)
    # alpha_3 = 48.0
    # beta_3 = 0.015
    # F_XSE = beta_3 * (exp(alpha_3 * A[:, 24]) - 1.0)

zero_dimentional()
# plotting_1_param()
# plotting_2_param()
