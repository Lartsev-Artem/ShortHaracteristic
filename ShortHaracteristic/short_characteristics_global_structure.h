#pragma once
#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "short_characteristics_headers.h"

typedef int IntId;
typedef double Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;

using namespace std;

struct Normals {
	std::vector<Vector3> n;
	Normals() {
	}

	Normals(const int size) {
		n.resize(size);
	}
};

struct cell {
	//int id;  - номер в массиве
	std::vector<Vector3> nodes_value;
	std::vector<int> neighbours_id_face;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3(-666, -666, -666));
		neighbours_id_face.resize(4, -5);
	}
};

#define PI 3.14159265358979323846
const double eps = 1e-10;

extern Vector3 start_point_plane_coord;   // начало координат плоскости
extern Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
extern Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

extern Matrix3	straight_face;  // 3 узла интерпол€ции
extern Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
extern vtkDataArray* density;
extern vtkDataArray* absorp_coef;
extern vtkDataArray* rad_en_loose_rate;

extern Type square_surface;  // площадь поверхности дискретной 

extern Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

extern int num_cur_direction; // номер текущего направлени€
extern Vector3 cur_direction;

extern int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности


extern const Vector3 center_point;
extern const Type inner_radius; // радиус внутренней сферы (с запасом)

#endif