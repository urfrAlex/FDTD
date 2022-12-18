import matplotlib.pyplot as plt
import numpy as np
import os
import xml.etree.ElementTree as et


if __name__== '__main__':
    x = np.arange(-100, 100, 1)
    y = -np.cos(x)*np.cos(np.pi)*np.exp(-(x-np.pi)**2)

    x_str = [str(j) for j in x]
    y_str = [str(j) for j in y]

    y_graph = -np.cos(x)*np.cos(np.pi)*np.exp(-(x-np.pi)**2)
    plt.plot(x, y_graph)
    plt.show()

    if not os.path.isdir("results"):
        os.mkdir("results")

    data = et.Element ('data')
    for i in range(len(x_str)):
        xdata = et.SubElement (data, 'xdata')
        x_res = et.SubElement(xdata, 'x').text = str(x_str[i])

    for i in range(len(x_str)):
        ydata = et.SubElement (data, 'ydata')
        y_res = et.SubElement(ydata, 'y').text = str(y_str[i])

    os.chdir(os.path.join(os.getcwd(), 'results'))

    xml_file = et.ElementTree(data)
    xml_file.write('results.xml')
