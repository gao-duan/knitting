#include <sstream>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <Windows.h>

#include "CImg.h"
#include "image.h"
#undef max
#undef min


const int nail_number = 315;
const int thread_path_length = 4000;
const int resolution = 200;
const Float unit_length = 0.05 / resolution; // 50mm / resolution
const Float Pi = std::acos(-1);
ColorMode color_mode = ColorMode::GRAY;
const int vis_resolution = 4096;
const int MAX_ITERS = 1000000;


int generate_random_nail() {
	return rand() % nail_number;
}


int generate_candidate_node() {
	return rand() % thread_path_length;
}

void LoadFromFile(std::vector<int>& thread_path, const std::string& name) {
	std::ifstream in(name);
	thread_path.clear();
	int x;
	while (in >> x) {
		thread_path.push_back(x);
	}

	for(int i = thread_path.size(); i < thread_path_length;++i) thread_path.push_back(generate_random_nail());
	in.close();
}

void DumpThread(const std::string& name, const std::vector<int>& thread_path) {
	std::ofstream out(name);
	for (auto t : thread_path) {
		out << t << std::endl;
	}
	out.close();
}


void init_path(std::vector<int>& thread_path) {
	thread_path.clear();
	for (int i = 0; i < thread_path_length; ++i) thread_path.push_back(generate_random_nail());
}

std::pair<int, int> get_coordinate_of_nail(int idx, int res = resolution) {
	Float angle = 2.0 * Pi / nail_number;

	Float x = std::cos(angle * idx) + 1;
	Float y = 1 - std::sin(angle * idx);

	x *= res / 2;
	y *= res / 2;

	int int_x = (int)x;
	int int_y = (int)y;
	if (int_x >= res) int_x = res - 1;
	if (int_y >= res) int_y = res - 1;
	

	return std::make_pair(int_x, int_y);
}

void update_current_img(const std::vector<int>& thread_path, Image& image) {
	image.Clear(0);

#pragma omp parallel for num_threads(8)
	for (int i = 0; i < thread_path_length - 1; ++i) {
		std::pair<int, int> p1 = get_coordinate_of_nail(thread_path[i]);
		std::pair<int, int> p2 = get_coordinate_of_nail(thread_path[i + 1]);
		image.DrawLineCount(p1, p2);
	}

	const int center = resolution / 2;

#pragma omp parallel for num_threads(8)
	for (int j = 0; j < resolution; ++j) {
		for (int i = 0; i < resolution; ++i) {
			Float dist = (i - center) * (i - center) + (j - center) * (j - center);
			if (dist > center * center - 2 * center * resolution * 0.005) {
				image.Get(i, j) = 0;
			}
		}
	}
	
	const Float min = 1.0 - image.Max();
	const Float _tmp = 1.0 / (1.0 - min);
	const Float max = image.Max();

#pragma omp parallel for num_threads(8)
	for (int j = 0; j < resolution; ++j) {
		for (int i = 0; i < resolution; ++i) {
			Float tmp = max - image.Get(i, j);
			image.Get(i, j) = tmp;
		}
	}
}

Float compute_error(const Image& target_image, const Image& image) {
	Float inv_mean1 = 1.0 / image.Sum(); 
	Float inv_mean2 = 1.0 / target_image.Sum(); 
	Float error = 0;
	
	const int center = resolution / 2;
#pragma omp parallel for reduction(+:error) num_threads(8)
	for (int j = 0; j < resolution; ++j) {
		for (int i = 0; i < resolution; ++i) {
			Float dist = (i - center) * (i - center) + (j - center) * (j - center);
			if (dist > center * center - 2 * center * resolution * 0.005) {
				continue;
			}
			Float tmp = image.Get(i, j) * inv_mean1  - target_image.Get(i, j) * inv_mean2;
			error += tmp *tmp;
		}
	}
	return error;
}

void Visualize(const std::string& dst_name, const std::vector<int>& thread_path) {
	Image vis_buffer(vis_resolution, vis_resolution);

	vis_buffer.Clear(1);

#pragma omp parallel for num_threads(8)
	for (int i = 0; i < thread_path_length - 1; ++i) {
		std::pair<int, int> p1 = get_coordinate_of_nail(thread_path[i], vis_resolution);
		std::pair<int, int> p2 = get_coordinate_of_nail(thread_path[i + 1], vis_resolution);
		vis_buffer.DrawLineBinary(p1, p2);
	}

	const int center = vis_resolution / 2;

#pragma omp parallel for num_threads(8)
	for (int j = 0; j < vis_resolution; ++j) {
		for (int i = 0; i < vis_resolution; ++i) {
			Float dist = (i - center) * (i - center) + (j - center) * (j - center);
			if (dist > center * center - 2 * center * vis_resolution * 0.005) {
				vis_buffer.Get(i, j) = 0;
			}
		}
	}
	vis_buffer.Write(dst_name, resolution, resolution);

}

void random_permute(const std::string& target_name, const std::string& dst_path, bool restart = false, const std::string& checkpoint_name="") {
	std::vector<int> thread_path;

	if (restart || checkpoint_name.empty()) {
		std::cout << "Random initial thread path..." << std::endl;
		init_path(thread_path);
	}
	else {
		std::cout << "Load thread path from " << checkpoint_name << std::endl;

		LoadFromFile(thread_path, checkpoint_name);
	}

	std::cout << "Load target image " << target_name << std::endl;
	Image target_image(resolution, resolution);
	target_image.Load(target_name, color_mode);
	
	std::ofstream log(dst_path + "log.txt");
	Float current_error = std::numeric_limits<Float>::infinity();
	int mutation = 0;

	std::cout << "Running optimizing..." << std::endl;

	Image image_buffer(resolution, resolution);  // update_current_img

	
	update_current_img(thread_path, image_buffer);	// initial image for first time 
	Image current_image(resolution, resolution);

	int iter = 0;
	do {
		// 1. randomly pick one candidate node and nail
		//    change node => candidate_nail
		const int node = generate_candidate_node();
		int current_nail = thread_path[node];

		int candidate_nail = generate_random_nail();
		thread_path[node] = candidate_nail;

		current_image = image_buffer;					// backup current image
		update_current_img(thread_path, image_buffer);	   // try to update

		Float error = compute_error(target_image, image_buffer);

		// apply mutation
		if (error < current_error) {
			current_error = error;
			mutation++;
		}
		// recover 
		else {
			thread_path[node] = current_nail;
			image_buffer = current_image;
		}		

		// print log and visualize
		if (iter % 10000 == 0) {
			//image.Write(dst_path + std::to_string(iter) + ".bmp");
			DumpThread(dst_path + std::to_string(iter) + ".txt", thread_path);
			std::cout << "iter: " << iter << " ; error: " << current_error <<"; mutation: " << mutation  << std::endl;
			log << "iter: " << iter << " ; error: " << current_error << "; mutation: " << mutation << std::endl;
		}		
		if (iter % 100000 == 0 || (iter % 1000 == 0 && iter < 10000) || (iter % 10000 == 0 && iter < 100000)) {
			Visualize(dst_path + "vis_" + std::to_string(iter) + ".png", thread_path);
		}
		iter += 1;

	} while (iter < MAX_ITERS);

	DumpThread(std::to_string(iter) + ".txt", thread_path);
	std::cout << "iter: " << iter << " ; error: " << current_error << "; mutation: " << mutation << std::endl;
	log << "iter: " << iter << " ; error: " << current_error << "; mutation: " << mutation << std::endl;
	Visualize(dst_path + "vis_final.png", thread_path);

	log.close();
}

void load_visualize(const std::string& checkpoint_name, const std::string& dst_name) {
	std::cout << "Load thread path from " << checkpoint_name << std::endl;
	std::vector<int> thread_path;
	LoadFromFile(thread_path, checkpoint_name);
	Visualize(dst_name, thread_path);
}

void CreateFolder(const char* path)
{
	if (!CreateDirectory(path, NULL))
	{
		return;
	}
}

int main(int argc, char** argv) {
	// input guide image
	const std::string target_name = "../dog.jpg";

	// output folder 
	const std::string dst_path = "../results/dog_gray2/";
	CreateFolder(dst_path.c_str());


	bool restart = true; // `false` to restore from previous checkpoint
	const std::string checkpoint = "../results/dog_gray/0.txt";

	random_permute(target_name, dst_path, restart, checkpoint);
	return 0;
}