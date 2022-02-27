#include "short_characteristics_main.h"
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
vtkDataArray* density;
vtkDataArray* absorp_coef;
vtkDataArray* rad_en_loose_rate;

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int num_cur_direction; // номер текущего направлени€
Vector3 cur_direction;

int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности

const Vector3 center_point(0, 0, 0);
const Type inner_radius = 0.51; // радиус внутренней сферы (с запасом)

int main(int argc, char* argv[])
{
	std::string name_file_settings = "";
	int max_number_of_iter = 1;
	Type accuracy = 1e-2;

	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file.txt";
	else
		name_file_settings = argv[1];
	if (argc > 2)
		max_number_of_iter = std::stoi(argv[2]);

	cout << "Max_number_of_iter= " << max_number_of_iter << '\n';

	size_t class_file_vtk;
	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string name_file_graph;
	std::string out_file_grid_vtk;
	std::string out_file_E1d;
	std::string name_file_normals;

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, name_file_graph, out_file_E1d,
		name_file_normals)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}


	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();  


	Type _clock = -omp_get_wtime();
	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";


	vector<Vector3> directions;
	vector<Type> squares;
	
	_clock = -omp_get_wtime();
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface);
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

	std::vector<cell> nodes_value;
	{
		std::vector<int> all_pairs_face;
		_clock = -omp_get_wtime();
		//int count_unique_face = FindNeighborsPairFace(unstructured_grid, all_pairs_face);
		FindNeighborsPairFaceAndBoundaries(unstructured_grid, all_pairs_face);

		InitNodesValue(all_pairs_face, nodes_value);

		_clock += omp_get_wtime();
	}
	std::cout << "\n Finding time of the all_pairs_face: " << _clock << "\n";


	const int count_directions = directions.size();
	const int count_cells = unstructured_grid->GetNumberOfCells();

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

	// ”пор€доченные индексы €чеек по данному направлению
	vector<IntId> sorted_id_cell(count_cells);

	Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/

	std::vector<Type> Illum(count_cells * count_directions, 0);
	std::vector<Type> Illum2(count_cells * count_directions, 0);


	std::vector<Normals> normals;
	ReadNormalFile(name_file_normals, normals);

	std::vector<Vector3> centers;
	FindAllCenterOfTetra(unstructured_grid, centers);

	int count = 0;
	Type norm = 0;
	Vector3 direction;

	ofstream ofile;
	ofile.open("File_with_Logs.txt");
	
	do {
		Type _clock = -omp_get_wtime();
		{

			/*---------------------------------- далее FOR по направлени€м----------------------------------*/
			for (int num_direction = 0; num_direction < count_directions; ++num_direction)
			{
				num_cur_direction = num_direction;
				direction = directions[num_direction];
				cur_direction = direction;

				ReadGraph(name_file_graph + to_string(num_direction) + ".txt", sorted_id_cell);
				ResetNodesValue(nodes_value);

				Vector3 x;
				Vector3 x0;
				int num_cell;
				const int count_cells = unstructured_grid->GetNumberOfCells();

				int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   

				/*---------------------------------- далее FOR по €чейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h) {
					num_cell = sorted_id_cell[h];					

					SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);
					FindInAndOutFaces(direction, num_cell, normals, face_state);

					for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
						if (!face_state[num_out_face])  // выход€щие грани
							GetNodes(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), num_out_face, vertex_tetra, face_state, direction,
								nodes_value,  Illum2, directions, squares);
					}

					Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), centers[num_cell], direction, vertex_tetra,
						nodes_value, Illum2, directions, squares);

					Illum[num_direction * count_cells + num_cell] = I_k_dir;
				}
				/*---------------------------------- конец FOR по €чейкам----------------------------------*/

				//std::cout << "End direction number: " << num_direction << '\n';
			}
			/*---------------------------------- конец FOR по направлени€м----------------------------------*/
		}
		Illum.swap(Illum2);
		
		_clock += omp_get_wtime();
		count++;
		norm = NormIllum(Illum, Illum2);
		std::cout << "Error:= " << norm << '\n';
		std::cout << "Time of iter: " << _clock << '\n';
		std::cout << "End iter_count number: " << count << '\n';

		ofile << "Error:= " << norm << '\n';
		ofile << "Time of iter: " << _clock << '\n';
		ofile << "End iter_count number: " << count << '\n';

	} while (norm > accuracy && count < max_number_of_iter);


	Illum.swap(Illum2);
	ofile.close();

	vector<Type> energy(count_cells);
	MakeEnergy(Illum, squares, square_surface, energy);
	std::cout << "energy\n";

	WriteFileSolution(out_file_grid_vtk, Illum, energy, unstructured_grid);
	return 0;
}


