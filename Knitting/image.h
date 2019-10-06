#pragma once

#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <algorithm>
#include <omp.h>
#include "CImg.h"
#undef max
#undef min

typedef double Float;

enum class ColorMode {GRAY, LUMIANCE, R, G, B, AVG};

struct Image {
	// Image() {}

	Image(int w, int h, Float default_val = 0):W(w), H(h) {
		data = new Float* [H];
		for (int i = 0; i < H; ++i) data[i] = new Float[W];
		Clear(default_val);
	}

	Image(const Image& img) {
		for (int i = 0; i < H; ++i) 
			delete[] data[i];
		delete[] data;
		
		H = img.H;
		W = img.W;

		data = new Float * [H];
		for (int i = 0; i < H; ++i) data[i] = new Float[W];
#pragma omp parallel for num_threads(8)
		for (int i = 0; i < H; ++i) {
			for (int j = 0; j < W; ++j) {
				data[i][j] = img.data[i][j];
			}
		}

	}

	~Image() {
		for (int i = 0; i < H; ++i) delete[] data[i];
		delete[] data;
	}

	Float Get(int x, int y) const {
		return data[y][x];
	}
	Float& Get(int x, int y) {
		return data[y][x];
	}


	Image& operator=(const Image& img) {
		if (H != img.H || W != img.W) {
			for (int i = 0; i < H; ++i)
				delete[] data[i];
			delete[] data;

			H = img.H;
			W = img.W;

			data = new Float * [H];
			for (int i = 0; i < H; ++i) data[i] = new Float[W];
		}
	
#pragma omp parallel for num_threads(8)
		for (int i = 0; i < H; ++i) {
			for (int j = 0; j < W; ++j) {
				data[i][j] = img.data[i][j];
			}
		}
		return *this;
	}

	void Clear(Float default_val) {
#pragma omp parallel for num_threads(8)
		for (int i = 0; i < H; ++i) {
			for (int j = 0; j < W; ++j) {
				data[i][j] = default_val;
			}
		}
	}

	void UpdatePixel(int x, int y, Float length) {
#pragma omp atomic 
		data[y][x] += length;
	}

	void DrawLineLength(const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
		int x1 = p1.first, y1 = p1.second;
		int x2 = p2.first, y2 = p2.second;

		if (x1 == x2) {
			int y_min = std::min(y1, y2);
			int y_max = std::max(y1, y2);
			for (int y = y_min; y <= y_max; ++y) {
				UpdatePixel(x1, y, 1);
			}
			return;
		}
		if (y1 == y2) {
			int x_min = std::min(x1, x2);
			int x_max = std::max(x1, x2);
			for (int x = x_min; x <= x_max; ++x) {
				UpdatePixel(x, y1, 1);

			}
			return;
		}

		Float k = (y2 - y1) / (x2 - x1);

		if ((x1 <= x2 && y1 <= y2) || (x1 >= x2 && y1 >= y2)) {
			int x_min = std::min(x1, x2);
			int x_max = std::max(x1, x2);
			int y_min = std::min(y1, y2);
			int y_max = std::max(y1, y2);

			Float ky = y_max - y_min + 1;
			Float kx = x_max - x_min + 1;

			int x = x_min, y = y_min;
			Float xx = x, yy = y;

			while (x <= x_max && y <= y_max) {
				int tmp_x = x, tmp_y = y;
				Float delta_x = x + 1 - xx;
				Float delta_y = y + 1 - yy;
				Float tx = delta_x / kx;
				Float ty = delta_y / ky;

				Float len;
				if (tx - ty > 1e-5) { // UP
					y += 1;
					len = std::sqrt((y - yy) * (y - yy) + (ty * kx) * (ty * kx));

					xx += ty * kx;
					yy = y;
				}
				else if (ty - tx > 1e-5) { // RIGHT
					x += 1;
					len = std::sqrt((x - xx) * (x - xx) + (tx * ky) * (tx * ky));

					xx = x;
					yy += tx * ky;

				}
				else { // UP-RIGHT
					x += 1;
					y += 1;
					len = std::sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));

					xx = x;
					yy = y;

				}


				UpdatePixel(tmp_x, tmp_y, len);

			}
		}
		else {
			// [xmin, ymax] => [xmax, ymin]
			int x_min = std::min(x1, x2);
			int x_max = std::max(x1, x2);
			int y_min = std::min(y1, y2);
			int y_max = std::max(y1, y2);

			int x = x_min, y = y_max;
			Float xx = x, yy = y;

			Float ky = y_max - y_min + 1;
			Float kx = x_max - x_min + 1;

			while (x <= x_max && y >= y_min) {
				int tmp_x = x, tmp_y = y;
				Float delta_x = x + 1 - xx;
				Float delta_y = yy - (y - 1);
				Float tx = delta_x / kx;
				Float ty = delta_y / ky;

				Float len;
				if (tx - ty > 1e-5) { // DOWN
					y -= 1;
					len = std::sqrt((y - yy) * (y - yy) + (ty * kx) * (ty * kx));

					xx += ty * kx;
					yy = y;
				}
				else if (ty - tx > 1e-5) { // RIGHT
					x += 1;
					len = std::sqrt((x - xx) * (x - xx) + (tx * ky) * (tx * ky));

					xx = x;
					yy -= tx * ky;

				}
				else { // DOWN-RIGHT
					x += 1;
					y -= 1;
					len = std::sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));

					xx = x;
					yy = y;

				}


				UpdatePixel(tmp_x, tmp_y, len);

			}


		}



	}

	void DrawLineCount(const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
		int x1 = p1.first, y1 = p1.second;
		int x2 = p2.first, y2 = p2.second;

		int dx, dy, i, e;
		int incx, incy, inc1, inc2;
		int x, y;

		dx = x2 - x1;
		dy = y2 - y1;

		if (dx < 0) dx = -dx;
		if (dy < 0) dy = -dy;
		incx = 1;
		if (x2 < x1) incx = -1;
		incy = 1;
		if (y2 < y1) incy = -1;
		x = x1; y = y1;
		if (dx > dy) {
			UpdatePixel(x, y, 1);
			e = 2 * dy - dx;
			inc1 = 2 * (dy - dx);
			inc2 = 2 * dy;
			for (i = 0; i < dx; i++) {
				if (e >= 0) {
					y += incy;
					e += inc1;
				}
				else
					e += inc2;
				x += incx;
				UpdatePixel(x, y, 1);
			}

		}
		else {
			UpdatePixel(x, y, 1);

			e = 2 * dx - dy;
			inc1 = 2 * (dx - dy);
			inc2 = 2 * dx;
			for (i = 0; i < dy; i++) {
				if (e >= 0) {
					x += incx;
					e += inc1;
				}
				else
					e += inc2;
				y += incy;
				UpdatePixel(x, y, 1);
			}
		}

	}

	void DrawLineBinary(const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
		int x1 = p1.first, y1 = p1.second;
		int x2 = p2.first, y2 = p2.second;

		int dx, dy, i, e;
		int incx, incy, inc1, inc2;
		int x, y;

		dx = x2 - x1;
		dy = y2 - y1;

		if (dx < 0) dx = -dx;
		if (dy < 0) dy = -dy;
		incx = 1;
		if (x2 < x1) incx = -1;
		incy = 1;
		if (y2 < y1) incy = -1;
		x = x1; y = y1;
		if (dx > dy) {
			data[y][x] = 0;
			e = 2 * dy - dx;
			inc1 = 2 * (dy - dx);
			inc2 = 2 * dy;
			for (i = 0; i < dx; i++) {
				if (e >= 0) {
					y += incy;
					e += inc1;
				}
				else
					e += inc2;
				x += incx;
				data[y][x] = 0;
			}

		}
		else {
			//data[y][x] += 1.0 / 255;
			data[y][x] = 0;
			e = 2 * dx - dy;
			inc1 = 2 * (dx - dy);
			inc2 = 2 * dx;
			for (i = 0; i < dy; i++) {
				if (e >= 0) {
					x += incx;
					e += inc1;
				}
				else
					e += inc2;
				y += incy;
				//data[y][x] += 1.0 / 255;
				data[y][x] = 0;
			}
		}

	}

	Float Max() const {
		Float shared_max = -1;
#pragma omp parallel 
		{
			Float max_value = -1;
#pragma omp for nowait

			for (int i = 0; i < H; ++i) {
				for (int j = 0; j < W; ++j) {
					max_value = std::max(max_value, data[i][j]);
				}
			}

#pragma omp critical 
			{
				shared_max = std::max(shared_max, max_value);
			}
		}
		return shared_max;
	}
	Float Min() const {
		Float shared_min = std::numeric_limits<Float>::max();
#pragma omp parallel 
		{
			Float min_value = std::numeric_limits<Float>::max();
#pragma omp for nowait

			for (int i = 0; i < H; ++i) {
				for (int j = 0; j < W; ++j) {
					min_value = std::min(min_value, data[i][j]);
				}
			}

#pragma omp critical 
			{
				shared_min = std::min(shared_min, min_value);
			}
		}
		return shared_min;
	}


	Float Mean() const {
		return Sum() / Float(W * H);
	}
	Float Sum() const {
		Float sum = 0;
#pragma omp parallel for reduction(+:sum) num_threads(8)
		for (int i = 0; i < H; ++i) {
			for (int j = 0; j < W; ++j) {
				sum += data[i][j];
			}
		}
		return sum;
	}

	void Write(const std::string& name, int _W=-1, int _H=-1) const {
		cimg_library::CImg<Float> image(W, H, 1, 1);
		for (int j = 0; j < H; j++) {
			for (int i = 0; i < W; i++) {
				Float c = data[j][i];
				c = c > 0 ? c : 0;
				c = c < 1 ? c : 1;
				c = c * 255;

				image(i, j, 0, 0) = c;

			}
		}
		if (_W > 0 && _H > 0) {
			image = image.resize(_W, _H, -100, -100, 6);
		}

		image.save(name.c_str());
	}
	void Load(const std::string& name, ColorMode mode = ColorMode::LUMIANCE) const {
		cimg_library::CImg<Float> image(name.c_str());
		if (image.width() != W || image.height() != H) {
			printf("Resize input image from [%d, %d] to [%d, %d].\n", image.width(), image.height(), W, H);
			image = image.resize(W, H, -100, -100, 6);
		}

		for (int j = 0; j < H; j++) {
			for (int i = 0; i < W; i++) {
				Float r = image(i, j, 0, 0) / 255.0;
				Float g = image(i, j, 0, 1) / 255.0;
				Float b = image(i, j, 0, 2) / 255.0;

				Float pixel_value;
				switch (mode)
				{
					case ColorMode::GRAY:
						pixel_value = 0.2989 * r + 0.5870 * g + 0.1140 * b;
						break;
					case ColorMode::LUMIANCE:
						pixel_value = 0.212671 * r + 0.715160 * g + 0.072169 * b;
						break;
					case ColorMode::R:
						pixel_value = r;
						break;
					case ColorMode::G:
						pixel_value = g;
						break;
					case ColorMode::B:
						pixel_value = b;
						break;
					case ColorMode::AVG:
						pixel_value = (r + g + b) / 3.0;
						break;
					default:
						pixel_value = 0.212671 * r + 0.715160 * g + 0.072169 * b;
						break;
				}
				data[j][i] = pixel_value;
			}
		}
	}

	Float** data;
	int W, H;
};