import numpy as np
from prodes.calculations.standard_equations import distance


class Grid():
    """An grid object which is used to divide 3D data into gridcells"""

    def __init__(self, size):
        self.size = size
        self.__cells = None
        

    @property
    def cells(self):
        """Calls the cells and checks if they are already constructed"""

        if self.__cells is None:
            raise NameError("No cells have been constructed")

        else:  
            return self.__cells

    def construct_cells(self, content):
        """constructs the cells and produces cells"""
        from prodes.core.point import Cell

        shape = self.scout_box_geometries(content)
        center = self.find_center(content)
        x_start, y_start, z_start = self.find_start(center, shape)
        x_cells, y_cells, z_cells = shape
    
        grid = []
        for z_number in range(z_cells):
            z = z_start + self.size * z_number            
            plane = []
            for y_number in range(y_cells):
                line = []
                y = y_start + self.size * y_number  
                for x_number in range(x_cells):
                    x = x_start + self.size * x_number
                    cell = Cell(self.size, x, y, z)
                    line.append(cell)
                plane.append(line)
            grid.append(plane)
        self.__cells = np.array(grid)

    def fill_cells(self, to_insert):
        """fills the grid cells with the components"""

        grid = self.cells
        for point in to_insert:
            x, y, z = self.in_which_cell(point)
            cell = grid[z][y][x]
            cell.add_content(point)

    def scout_box_geometries(self, content, factor=1):
        """scouts the different axis lengths of the grid"""
        
        shape = []
        for coordinate in ["x", "y", "z"]:
            extremes = self.find_extremes(np.array([getattr(component, coordinate) for component in content])) * factor
            shape.append(self.required_cells(*extremes))
        
        return shape

    def find_extremes(self, array):
        """finds the extremes of a list of numbers"""
        
        return np.amin(array), np.amax(array)

    def find_center(self, content):
        """finds the center of all components in the grid"""

        center = []
        for coordinate in ["x", "y", "z"]:
            extreme_min, extreme_max = self.find_extremes(np.array([getattr(component, coordinate) for component in content]))
            center.append((extreme_max + extreme_min)/ 2)
        
        return center

    def required_cells(self, minimum, maximum):
        """Calculates the number of required cells to based on the cell size"""
        
        from math import ceil, floor
        minimum, maximum = floor(minimum), ceil(maximum)
        return ceil((maximum- minimum)/self.size)

    def find_start(self, center, shape):
        """finds the x, y and z coordinates of the starting cells"""

        starts = []
        for center_coordiante, n_cells in zip(center, shape):
            start = center_coordiante - n_cells / 2 * self.size + self.size / 2
            starts.append(start)
        return starts

        
    def in_which_cell(self, point):
        """Finds the coordinates in which cell a point should be present"""

        from math import floor

        corner_cell  = self.cells[0, 0, 0]
        which_cell = []
        for coordinate in ["x", "y", "z"]:
            point_coordinate = getattr(point, coordinate)
            origin_coordinate = getattr(corner_cell, coordinate) - self.size / 2
            which_cell.append(floor((point_coordinate - origin_coordinate) / self.size))

        return which_cell

    def find_surrounding_cells(self, cell):
        """Gets the gridcells surounding a specific cell"""

        grid = self.cells
        cell_x, cell_y, cell_z = self.in_which_cell(cell)
        z_lim, y_lim, x_lim = [numb - 1 for numb in grid.shape]
        sub_enviroment = []
        directions = [-1, 0, 1]
        for z in directions:
            to_add_z = cell_z+ z
            if to_add_z >= 0 and to_add_z <= z_lim:
                for y in directions:
                    to_add_y = cell_y + y
                    if to_add_y >= 0 and to_add_y <= y_lim:
                        for x in directions:
                            to_add_x = cell_x + x
                            if to_add_x >= 0 and to_add_x <= x_lim:
                                sub_enviroment.append(grid[to_add_z, to_add_y, to_add_x])
        return np.array(sub_enviroment)

    def grid_content(self, *content_filter, cells=None):
        """Returns the content of the cells
        takes strings representing object names as arguments to filter for specific objects
        optional argument cells allow for to get the content of a subset of cells"""

        if cells is None:
            cells = self.cells
            cells = cells.flatten()

        if len(content_filter) == 0:
            content = np.array([point for cell in cells for point in cell.content])
        else:
            content = np.array([point for cell in cells for point in cell.content if type(point).__name__ in content_filter])

        return content
        


def property_points_on_surface(grid, to_consider="Surface_point"):
    """Finds the property points on a surface"""

    from prodes.core.point import Property_point

    for cell in grid.cells.flatten():
        if len(cell.filtered_content(to_consider)) > 0:
            start = cell
            break

    surface_points = np.empty([0])
    done_cells = np.empty([0])
    todo_cells = np.array([start])

    while todo_cells.size > 0:
        next_todo = np.empty([0])
        for cell in todo_cells:
            surface_point = Property_point(cell.x, cell.y, cell.z)
            surface_points = np.concatenate([surface_points, np.array([surface_point])])
            done_cells = np.concatenate([done_cells, np.array([cell])])
            
            surrounding = grid.find_surrounding_cells(cell)

            for neighbour in surrounding:
                if neighbour.empty == cell.empty and neighbour not in next_todo and neighbour not in done_cells:
                    next_todo = np.concatenate([next_todo, np.array([neighbour])])

        todo_cells = np.array([cells for cells in next_todo if cells not in done_cells])

    return surface_points
