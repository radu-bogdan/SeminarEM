import meshio
import numpy as np

mesh = meshio.read(
    "Motor_Bosch_2d.vol",  # string, os.PathLike, or a buffer/open file
    file_format = "netgen",  # optional if filename is a path; inferred from extension
    # see meshio-convert -h for all possible formats
)

p = mesh.points
t = np.c_[mesh.cells_dict['triangle'],
          mesh.cell_data_dict['netgen:index']['triangle']]
e = np.c_[mesh.cells_dict['line'],
          mesh.cell_data_dict['netgen:index']['line']]