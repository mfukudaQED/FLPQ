import numpy as np
import sys
import os

BohrR = 0.529177249 # Angstrom


class Data_LPQ:
    """ class of data from md file """
    def __init__(self,
                 atomnum = None,
                 numbers = None,   # int
                 valence = None,   # int
                 positions_bohr = None,  # bohr, coordinate of atoms 
                 positions_ang = None,  # Ang, coordinate of atoms 
                 origin = None, # Ang
                 num_mesh = None,
                 mesh_vec = None, # Ang primitive vector for mesh
                 mesh_grid = None, # Ang
                 lattice = None, # Ang
                 dim_val = None,  # int
                 val = None,  # scalar, vector, tensor, etc.
                 small_mesh_grid = None, # Ang
                 small_val = None):

        pass

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

 
    def read_header(self,filename):
        with open(filename, 'r') as file:
    
            # Skip the header lines
            for _ in range(2):
                next(file)
            
            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            self.origin = [float(coord) for coord in line[1:4]]
            self.origin = [float(coord)*BohrR for coord in line[1:4]]
            
            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]

            self.lattice = []
            for i in range(3):
                self.lattice.append([])
                for j in range(3):
                    self.lattice[i].append(self.mesh_vec[i][j]*self.num_mesh[i])

    
            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            for _ in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])

            ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)
        return 
         
 
    def read_val(self,filename):

        with open(filename, 'r') as file:
            ### skip two comment lines ###
            line = next(file).strip().split()
            line = next(file).strip().split()

            line = next(file).strip().split()
            self.dim_val = len(line) - 3

        with open(filename, 'r') as file:
            ### skip two comment lines ###
            line = next(file).strip().split()
            line = next(file).strip().split()

            num_data = self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2]
            self.val = []
            for i in range(self.dim_val):
                self.val.append(np.empty(num_data).tolist())
    
            # Read the file using the next function until the end of the file
            idx = 0
            while True:
                try:
                    line = next(file).strip().split()
                    for value in line[3:]:
                        self.val[int(idx%(self.dim_val))][int(idx/(self.dim_val))] = float(value)
                        #print(self.val[idx])
                        idx += 1
    
                except StopIteration:
                    # End of file reached, break out of the loop
                    break

        return 

    def read_val_new(self,filename):

        with open(filename, 'r') as file:
            line = next(file).strip().split()
            self.dim_val = int(line[0])
            num_data = self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2]
            self.val = []
            for i in range(self.dim_val):
                self.val.append(np.empty(num_data).tolist())
    
            # Read the file using the next function until the end of the file
            idx = 0
            while True:
                try:
                    line = next(file).strip().split()
                    for value in line:
                        self.val[int(idx%self.dim_val)][int(idx/self.dim_val)] = float(value)
                        #print(self.val[idx])
                        idx += 1
    
                except StopIteration:
                    # End of file reached, break out of the loop
                    break

        return
      
    def set_mesh_grid(self):

        num_data = self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2]
        self.mesh_grid = np.empty((3,num_data))

        x = self.origin[0]
        y = self.origin[1]
        z = self.origin[2]
        idx = 0
        for i in range(self.num_mesh[0]):
            for j in range(self.num_mesh[1]):
                for k in range(self.num_mesh[2]):
                    for l in range(3):
                        self.mesh_grid[l][idx] = self.origin[l] + i*self.mesh_vec[0][l] +j*self.mesh_vec[1][l] +k*self.mesh_vec[2][l]
                    #if ( self.val[idx] > 1e-4 ):
                    #    print(self.mesh_grid[0][idx],self.mesh_grid[1][idx],self.mesh_grid[2][idx],self.val[idx])
                    idx += 1 
        return 

    
    def set_small_mesh_grid(self,points,mergin):

        points_max = []
        points_min = []
        for i in range(3):
            points_max.append(np.max(points[i]))
            points_min.append(np.min(points[i]))

        #mergin = 0.5 # Ang


        list_idx = np.zeros(len(self.mesh_grid[0]))

        for i in range(len(self.mesh_grid[0])):
            if (self.mesh_grid[0][i] < (points_max[0] + mergin)):
                if (self.mesh_grid[1][i] < (points_max[1] + mergin)):
                    if (self.mesh_grid[2][i] < (points_max[2] + mergin)):
                        if (self.mesh_grid[0][i] > (points_min[0] - mergin)):
                            if (self.mesh_grid[1][i] > (points_min[1] - mergin)):
                                if (self.mesh_grid[2][i] > (points_min[2] - mergin)):
                                    list_idx[i] = 1

        num_small_mesh_grid = int(np.sum(list_idx))
        self.small_mesh_grid = np.empty((3,num_small_mesh_grid))
        self.small_val = []
        for i in range(self.dim_val):
            self.small_val.append( np.empty(num_small_mesh_grid) )

        idx = 0
        for i in range(len(self.mesh_grid[0])):
            if (list_idx[i]==1):
                for j in range(3):
                    self.small_mesh_grid[j][idx] = self.mesh_grid[j][i]
                for k in range(self.dim_val):
                    self.small_val[k][idx] = self.val[k][i]
                idx += 1

        return 
    
    #def set_cell_index_old(self,points,mergin):
    #    
    #    lattice = np.array(self.lattice)
    #    origin = np.array(self.origin)
    #    transformation_matrix = np.linalg.inv(lattice.T)

    #    xmesh = len(points[0])
    #    ymesh = len(points[1])
    #    zmesh = len(points[2])
    #    num_mesh = xmesh*ymesh*zmesh

    #    cell_index = []
    #    for idx in range(num_mesh):
    #        cell_index.append([])
    #        for i in range(3):
    #            cell_index[idx].append(0)

    #    ndx = 0
    #    for idx in range(xmesh):
    #        for jdx in range(ymesh):
    #            for kdx in range(zmesh):
    #                point = np.array([points[0][idx],points[1][jdx],points[2][kdx]])
    #                transformed_point = np.dot(transformation_matrix, (point - origin))
    #                #print(transformed_point)

    #                for i,coord in enumerate(transformed_point):
    #                    flag = True
    #                    while (flag):
    #                        if coord < 0:
    #                             coord = coord + 1
    #                             cell_index[ndx][i] = int(cell_index[ndx][i])-1
    #                        elif coord >= 1:
    #                             coord = coord - 1
    #                             cell_index[ndx][i] = int(cell_index[ndx][i])+1
    #                        else:
    #                             flag = False
    #                ndx +=1


    #    points_max = []
    #    points_min = []
    #    for i in range(3):
    #        points_max.append(np.max(points[i]))
    #        points_min.append(np.min(points[i]))

    #    #mergin = 0.5 # Ang


    #    list_idx = np.zeros(len(self.mesh_grid[0]))

    #    for i in range(len(self.mesh_grid[0])):
    #        if (self.mesh_grid[0][i] < (points_max[0] + mergin)):
    #            if (self.mesh_grid[1][i] < (points_max[1] + mergin)):
    #                if (self.mesh_grid[2][i] < (points_max[2] + mergin)):
    #                    if (self.mesh_grid[0][i] > (points_min[0] - mergin)):
    #                        if (self.mesh_grid[1][i] > (points_min[1] - mergin)):
    #                            if (self.mesh_grid[2][i] > (points_min[2] - mergin)):
    #                                list_idx[i] = 1

    #    num_small_mesh_grid = int(np.sum(list_idx))
    #    self.small_mesh_grid = np.empty((3,num_small_mesh_grid))
    #    self.small_val = np.empty((self.dim_val,num_small_mesh_grid))

    #    idx = 0
    #    for i in range(len(self.mesh_grid[0])):
    #        if (list_idx[i]==1):
    #            for j in range(3):
    #                self.small_mesh_grid[j][idx] = self.mesh_grid[j][i]
    #            for k in range(self.dim_val):
    #                self.small_val[k][idx] = self.val[k][i]
    #            idx += 1

    #    return 
    
    def set_cell_index(self,points,points_meshgrid,mergin):
        
        lattice = np.array(self.lattice)
        origin = np.array(self.origin)
        transformation_matrix = np.linalg.inv(lattice.T)

        cell_index = []
        for idx in range(len(points)):
            cell_index.append([])
            for i in range(3):
                cell_index[idx].append(0)

        for idx in range(len(points)):
            point = np.array([points[idx][0],points[idx][1],points[idx][2]])
            transformed_point = np.dot(transformation_matrix, (point - origin))
            #print(transformed_point)

            for i,coord in enumerate(transformed_point):
                flag = True
                while (flag):
                    if coord < 0:
                         coord = coord + 1
                         cell_index[idx][i] = int(cell_index[idx][i])-1
                    elif coord >= 1:
                         coord = coord - 1
                         cell_index[idx][i] = int(cell_index[idx][i])+1
                    else:
                         flag = False

        points_max = []
        points_min = []
        for i in range(3):
            points_max.append(np.max(points_meshgrid[i]))
            points_min.append(np.min(points_meshgrid[i]))

        #mergin = 0.5 # Ang


        list_idx = np.zeros(len(self.mesh_grid[0]))

        for i in range(len(self.mesh_grid[0])):
            if (self.mesh_grid[0][i] < (points_max[0] + mergin)):
                if (self.mesh_grid[1][i] < (points_max[1] + mergin)):
                    if (self.mesh_grid[2][i] < (points_max[2] + mergin)):
                        if (self.mesh_grid[0][i] > (points_min[0] - mergin)):
                            if (self.mesh_grid[1][i] > (points_min[1] - mergin)):
                                if (self.mesh_grid[2][i] > (points_min[2] - mergin)):
                                    list_idx[i] = 1

        num_small_mesh_grid = int(np.sum(list_idx))
        self.small_mesh_grid = np.empty((3,num_small_mesh_grid))
        self.small_val = np.empty((self.dim_val,num_small_mesh_grid))

        idx = 0
        for i in range(len(self.mesh_grid[0])):
            if (list_idx[i]==1):
                for j in range(3):
                    self.small_mesh_grid[j][idx] = self.mesh_grid[j][i]
                for k in range(self.dim_val):
                    self.small_val[k][idx] = self.val[k][i]
                idx += 1

        return 

    def interpolate_data(self,points,mesh_grid,val):
        from scipy.interpolate import griddata

        # points : (x, y, z)
        # val : scalar

        interpolated_data = []

        #print(len(mesh_grid[0]))
        #print(len(val[0]))
        #print(points)

        #points_meshgrid = np.meshgrid(points[0], points[1], points[2], indexing='ij', sparse=False)
        #xi = np.column_stack((points_meshgrid[0].ravel(), points_meshgrid[1].ravel(), points_meshgrid[2].ravel()))

        #print(xi)

        input_grid = np.column_stack((mesh_grid[0].ravel(), mesh_grid[1].ravel(), mesh_grid[2].ravel()))
        #print(points_meshgrid)
        for i in range(len(val)):
            interpolated_data.append( griddata(input_grid, val[i], points, method='linear') ) # for 3D

        return interpolated_data

    #def mesh_grid_extended(self,mesh_grid,val):
    #    return 

    #def plot_value_on_contour_3d_mayavi(self,points,val,val_contours):
    #    from mayavi import mlab
    #    
    #    # 3次元座標とスカラー量のデータを生成する
    #    x, y, z = points[0],points[1],points[2]
    #    scalar_field = val[0]
    #    scalar_field = np.reshape(scalar_field, x.shape)
    #    print(x)
    #    print(scalar_field)
    #    
    #    # isosurfaceを作成する
    #    mlab.figure()
    #    #mlab.contour3d(points[0],points[1],points[2], val, contours=[0.5])
    #    mlab.contour3d(x, y, z, scalar_field, contours=val_contours)
    #    
    #    ## 座標データを取得する
    #    isosurface = mlab.gcf().children[0]
    #    coords = isosurface.mlab_source.points.to_array()
    #    #print(coords)

    #    surf = mlab.surf(coords)

    #    # スカラー量のカラープロットを設定する
    #    surf.module_manager.scalar_lut_manager.data_range = (np.min(scalar_field), np.max(scalar_field))
    #    surf.module_manager.scalar_lut_manager.scalar_data = scalar_field.flatten()
    #    
    #    # カラーバーを表示する
    #    mlab.colorbar()

    #    # プロットを表示する
    #    #mlab.show()
    #    mlab.savefig("contour3d.pdf")


    def plot_value_on_contour_3d(self,points_meshgrid,val,val_contours):
        import matplotlib.pyplot as plt
        import matplotlib.colors
        
        # 3次元座標とスカラー量のデータを生成する
        x, y, z = points_meshgrid[0],points_meshgrid[1],points_meshgrid[2]
        alpha = val[0]

        # Create the contour plot
        contour = plt.contour3d(x, y, z, alpha, levels=val_contours)

        # Extract the coordinates of the isosurface
        isosurface = contour.collections[0].get_paths()[0]
        coordinates = isosurface.vertices

        # カラーマップの作成
        my_cmap = plt.get_cmap('viridis')
        
        # alpha値を[0,1]の範囲に正規化
        norm=matplotlib.colors.Normalize(vmin=alpha.min(), vmax=alpha.max())
        colors = my_cmap(norm(alpha))
        
        # フィギュアの作成
        fig = plt.figure(figsize =(14, 9))
        ax = plt.axes(projection ='3d')
        
        # プロットの作成
        surf = ax.plot_surface(x, y, z, facecolors=colors, edgecolor ='none')
        
        # カラーバーの設定
        sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
        sm.set_array(alpha)
        fig.colorbar(sm, shrink=0.5, aspect=5)
        
        ax.set_title('Surface plot with colormap based on alpha')
        
        # プロットの表示
        #plt.show()
        plt.savefig("contour3d.pdf")


def set_isosurface_skimage(num_mesh,meshgrid,val,val_contours):
    from skimage import measure
    
    # 3次元座標とスカラー量のデータを生成する
    x, y, z = meshgrid[0],meshgrid[1],meshgrid[2]

    data_volume = np.empty((num_mesh[0],num_mesh[1],num_mesh[2]))

    idx = 0
    for i in range(num_mesh[0]):
        for j in range(num_mesh[1]):
            for k in range(num_mesh[2]):
                data_volume[i,j,k] = val[idx]
                idx +=1
    #print(val[0])
    #print(data_volume)

    #points = []
    #for idx in range(len(x)):
    #    points.append([x[idx],y[idx],z[idx]])

    #selected_points = [point for i, point in enumerate(points) if alpha[i] >= val_contours[0]]
    #print(selected_points)
 
    # Isosurfaceの生成
    vertices, faces, normals, values = measure.marching_cubes(data_volume, level=val_contours[0])
    
    # Isosurface上の座標を取得
    coords = vertices
    print(coords)

    return coords

def set_isosurface_mayavi(meshgrid,val,val_contours):
    from mayavi import mlab
    
    # 3次元座標とスカラー量のデータを生成する
    x, y, z = meshgrid[0],meshgrid[1],meshgrid[2]
    alpha = val
 
    # Create the contour plot
    mlab.contour3d(x, y, z, alpha, contours=val_contours)
    
    ## 座標データを取得する
    isosurface = mlab.gcf().children[0]
    coords = isosurface.mlab_source.points.to_array()

    return coords

##########################
if __name__ == '__main__':
    
    filepath = os.getcwd()
    current_dir = os.path.basename(os.path.dirname(filepath))

    filename_cube_header = sys.argv[1]
    filename1 = sys.argv[2]
    filename2 = sys.argv[3]
    filename_frmsf = "data.frmsf"

    data_LPQ = Data_LPQ()
    data_LPQ.read_header(filename_cube_header)
    data_LPQ.read_val(filename1)

    data_LPQ2 = Data_LPQ()
    data_LPQ2.read_header(filename_cube_header)
    data_LPQ2.read_val(filename2)

    with open(filename_frmsf, 'w') as f:
        f.write(f"{data_LPQ.num_mesh[0]} {data_LPQ.num_mesh[1]} {data_LPQ.num_mesh[2]}")
        f.write("\n")
        f.write("1")
        f.write("\n")
        f.write("1") # nband
        f.write("\n")
        f.write(f"{data_LPQ.lattice[0][0]:.6e} {data_LPQ.lattice[0][1]:.6e} {data_LPQ.lattice[0][2]:.6e}")
        f.write("\n")
        f.write(f"{data_LPQ.lattice[1][0]:.6e} {data_LPQ.lattice[1][1]:.6e} {data_LPQ.lattice[1][2]:.6e}")
        f.write("\n")
        f.write(f"{data_LPQ.lattice[2][0]:.6e} {data_LPQ.lattice[2][1]:.6e} {data_LPQ.lattice[2][2]:.6e}")
        f.write("\n")
        idx = 0
        for i in range(data_LPQ.num_mesh[0]):
            for j in range(data_LPQ.num_mesh[1]):
                for k in range(data_LPQ.num_mesh[2]):
                    f.write(f"{data_LPQ.val[0][idx]:.3e}")
                    f.write("\n")
                    idx += 1
        idx = 0
        for i in range(data_LPQ.num_mesh[0]):
            for j in range(data_LPQ.num_mesh[1]):
                for k in range(data_LPQ.num_mesh[2]):
                    f.write(f"{data_LPQ2.val[0][idx]:.3e}")
                    f.write("\n")
                    idx += 1


    exit()

    data_LPQ.set_mesh_grid()

    val_contours = [0.00001]
    points = set_isosurface_skimage(data_LPQ.num_mesh, data_LPQ.mesh_grid, data_LPQ.val[0], val_contours)

    #print(points)
    exit()

    ### set points and mergin ###
    def set_arbitrary_points():
        x_values = np.linspace(-5, 5, 11)
        y_values = np.linspace(-5, 5, 11)
        z_values = np.linspace(-5, 5, 11)
        points = []
        for x in range(x_values):
            for y in range(y_values):
                for z in range(z_values):
                    points.append([x, y, z])

        return points

    def convert_points_to_meshgrid(points):
        for idx in range(len(points)):
            points_meshgrid_x = points[idx][0]
            points_meshgrid_y = points[idx][1]
            points_meshgrid_z = points[idx][2]

        points_meshgrid = [points_meshgrid_x, points_meshgrid_y, points_meshgrid_z]

        return points_meshgrid

    def set_input_points():
        x_values = np.linspace(-5, 5, 11)
        y_values = np.linspace(-5, 5, 11)
        z_values = np.linspace(-5, 5, 11)
        input_points = (x_values, y_values, z_values)
        #input_meshgrid = np.meshgrid(x_values, y_values, z_values, indexing='ij', sparse=True)

        return input_points

    def set_input_meshgrid():
        x_values = np.linspace(-5, 5, 11)
        y_values = np.linspace(-5, 5, 11)
        z_values = np.linspace(-5, 5, 11)
        input_meshgrid = np.meshgrid(x_values, y_values, z_values, indexing='ij', sparse=False)

        return input_meshgrid

    points = set_arbitrary_points()
    points_meshgrid = convert_points_to_meshgrid(points)
    #print(points)

    mergin = 0.5

    data_LPQ.set_cell_index(points,mergin)
    data_LPQ.set_small_mesh_grid(points,mergin)
    val = data_LPQ.interpolate_data(points,data_LPQ.small_mesh_grid,data_LPQ.small_val)
    val_contours = [0.00001]
    data_LPQ.plot_value_on_contour_3d(points,val[0],val_contours)
