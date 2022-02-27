#include "short_characteristics_logic_function.h"

int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани				
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}// x->координата узла на выходящей грани		}
	//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}

	int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];
	if(neighbor_id_face>=0)
	nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
		nodes_value[num_cur_cell].nodes_value[num_cur_out_face];
	return 0;
}

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0)) {

			Type s = (x - x0).norm();

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value);

			Type I = CurGetIllum(num_cur_cell, x0, s, I_x0, direction, illum_old, directions, squares);

			nodes_value[num_cur_cell].nodes_value[num_cur_out_face][num_node] = I;

			break;
		}

	}//for num_in_face

	return 0;
}

size_t Intersection(const std::vector<Vector3>& face, const Vector3& X0, const Vector3& n, std::vector<Type>& res) {
	/*Type* A = face[0];
	Type* B = face[1];
	Type* C = face[2];*/

	Vector3 A = face[0];
	Vector3 B = face[1];
	Vector3 C = face[2];

	Type a, b, c, d;
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	//cout << setprecision(16)<<"a= " << a << "\nb= " << b << "\nc= " << c << "\nd= " << d << '\n';

	t = -(a * X0[0] + b * X0[1] + c * X0[2] + d) / (a * n[0] + b * n[1] + c * n[2]);

	for (size_t i = 0; i < 3; i++)
		res[i] = (n[i] * t + X0[i]);

	return 0;
}

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value) {
	Type I_x0 = 0;
	
	if (nodes_value[num_cell].neighbours_id_face[num_in_face] == -1) {
		/*Граничные условия*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;
	}
	else if (nodes_value[num_cell].neighbours_id_face[num_in_face] == -2) {

		// внутренняя граница
		const Type Rsphere = 0.1;
		const Type R1disk = 0.1; 
		const Type R2disk = 0.2;

		std::vector<Vector3> curface(3);		

		std::vector<Type> res(3,0);
		
		//пересечение с диском
		{
			// точки задающие плоскость диска
			curface[0][0] = 1;
			curface[0][1] = 0;
			curface[0][2] = 0;

			curface[1][0] = 0;//0;
			curface[1][1] = 0.9928768384869221;//
			curface[1][2] = 0.11914522061843064;//;

			curface[2][0] = 2;//;
			curface[2][1] = 0;
			curface[2][2] = 0;// ;   // Wolfram

			// пересечние луча с плоскостью диска
			Intersection(curface, x, cur_direction, res);

			Type LocRes[3] = { 0,0,0 }; // точка пересечения в локальных координатах плоскости

			Type v1[3] = { 1, 0, 0 };
			Type v2[3] = { 0, -0.992877, -0.119145 }; // Wolfram
			//в плоскости
			for (size_t k = 0; k < 3; k++) {
				LocRes[0] += (res[k]) * v1[k];
				LocRes[1] += (res[k]) * v2[k];
			}

			Type A = pow(cur_direction[0], 2) + pow(cur_direction[1], 2) + pow(cur_direction[2], 2);
			Vector3 a = cur_direction;

			Type radical = 4 * pow((a[0] * (-1 + x[0]) + a[1] * x[1] + a[2] * x[2]), 2) -
				4 * A * (1 - 2 * x[0] + pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2) - pow(Rsphere, 2));

			if (radical >= 0) {// есть пересечение со сферой
				Type t = (a[0] - a[0] * x[0] - a[1] * x[1] - a[2] * x[2] -
					sqrt(radical) / 2) / A;

				Type inSphere[3] = { a[0] * t + x[0], a[1] * t + x[1], a[2] * t + x[2] };

				// не пересекает плоскость
				if (pow(LocRes[0] - center_point[0], 2) + pow(LocRes[1], 2) <= pow(R1disk, 2) || 
					pow(LocRes[0] - center_point[0], 2) + pow(LocRes[1], 2) >= pow(R2disk, 2))
					return 50; // 2;

				// с чем луч встречается раньше?
				Type LenSpehere = pow(x[0] - inSphere[0], 2) + pow(x[1] - inSphere[1], 2) + pow(x[2] - inSphere[2], 2);
				Type LenPlane = pow(x[0] - res[0], 2) + pow(x[1] - res[1], 2) + pow(x[2] - res[2], 2);

				if (LenSpehere > LenPlane)
					return 20; // 1;
				else
					return 50; // 2;
			}
			else if (pow(LocRes[0] - center_point[0], 2) + pow(LocRes[1], 2) < pow(R2disk, 2) &&
				(pow(LocRes[0] - center_point[0], 2) + pow(LocRes[1], 2) > pow(R1disk, 2)))
				return 20; // 1;

		}

		return 10;// I_x0;

	}
	else {

		if (nodes_value[num_cell].nodes_value[num_in_face][0] < -600)
			cout << "Num_dir: " << num_cur_direction << " CalculateIllumeOnInnerFace:  Undefine cell / in face:" << num_cell << " / " << num_in_face << " !!!\n";

		Vector3 x0_local;

		FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);
		Vector3 coef;// = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);

		switch (num_in_face) {
		case 3:
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);

			coef = GetInterpolationCoef(inclined_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = local_plane_x0[0] * coef[0] + local_plane_x0[1] * coef[1] + coef[2];

			/*I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];
			Vector3 coef = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);*/
			break;
		case 1:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 2:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 0:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

			break;
		}
		if (I_x0 < 0) {
			count_negative_interpolation++;
			return 0;
		}

		return I_x0;
	}
}

Type CurGetIllum(const int cur_id, const Vector3 x, const Type s, const Type I_node_prev, const Vector3& cur_direction,
	const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {
	// без интеграла рассеивания
		{
				Type Ie = 10;
				Type k = 10;
				if (x.norm() > 0.25) { Ie = 0; k = 1; }

				Type I;
				if (s > 1e-10)
					I = Ie * (1 - exp(-s * k)) + I_node_prev * exp(-s * k);
				else
					I = Ie * (1 + s * k) - I_node_prev * s * k;

				if (I < 0)
					I = 0;
				return I;
		}


		Type S = GetS(cur_id, cur_direction, illum_old, directions, squares);
		Type Ie = 10.;
		Type alpha = 5.;
		Type betta = 5.;
		Type k = alpha + betta;
		if (x.norm() > 0.3) {
			Ie = 0;
			alpha = 0.5;
			betta = 0.5;
			k = alpha + betta;
		}

		Type I = exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Ie * alpha + S * betta));
		I /= k;

		if (I < 0)
			I = 0;
		return I;
}

Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {
	/*Все грани должно быть определены*/
	Type value = -666;
	Vector3 x0;

	for (size_t i = 0; i < 4; i++) {

		IntersectionWithPlane(cur_cell->GetFace(i), center, direction, x0);
		if (InTriangle(num_cell, unstructuredgrid, cur_cell, i, x0)) {
			if ((center - x0).dot(direction) <= 0) continue;
			
			Type s = (center - x0).norm();
			Type I_x0 = CalculateIllumeOnInnerFace(num_cell, i, vertex_tetra, center, x0, nodes_value);

			value = CurGetIllum(num_cell, x0, s, I_x0, direction, illum_old, directions, squares);
			break;
		}
	}
	if (value < 0)
		return 0;
	return value;
}
