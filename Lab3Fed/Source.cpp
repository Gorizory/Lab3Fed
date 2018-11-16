#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <ctime>
#include <vector>
#include <string>

#define N 3
#define NU 0.4
#define N_ITER 3

using namespace std;

class Point {
public:
	double x, y, d;

	Point(double _x, double _y, double _d) {
		x = _x;
		y = _y;
		d = _d;
	}
};

vector<double> normalizeVector(vector<double> vec) {
	double inSum = 0.0;

	for (double j : vec) {
		inSum += pow(j, 2);
	}

	double norm = sqrt(inSum);

	for (unsigned j = 0; j < vec.size(); j++) {
		vec[j] = vec[j] / norm;
	}

	return vec;
}

double randDouble() {
	return (double) (rand() / (double) RAND_MAX) - 1.0;
}

vector<Point> createEllipce(double x0, double y0, double angle, double a, double b, double d) {
	double x, y, ellipseArea;
	vector<Point> vec;

	double angleRad = M_PI / 180.0 * angle;

	for (x = -50.0; x < 50.0; x += 0.5) {
		for (y = -50.0; y < 50.0; y += 0.5) {
			double _x = x + randDouble();
			double _y = y + randDouble();

			ellipseArea =
				pow((_x - x0) * cos(angleRad) + (_y - y0) * sin(angleRad), 2) / pow(a, 2) +
				pow((-(_x - x0) * sin(angleRad)) + (_y - y0) * cos(angleRad), 2) / pow(b, 2);

			if (ellipseArea <= 1.0) {
				vec.push_back(Point(_x, _y, d));
			}
		}
	}

	return vec;
}

class Grossberg {
private:
	vector<double> w;
public:
	Grossberg() {
		for (int j = 0; j < N; j++) {
			w.push_back(0.0);
		}
	}

	void normW() {
		w = normalizeVector(w);
	}

	vector<double> getW() {
		return w;
	}

	vector<double> teach(Point point) {
		vector<double> vec = { 0.0 };
		vector<double> delta = { 0.0 };

		vec.push_back(point.x);
		vec.push_back(point.y);

		vec = normalizeVector(vec);

		for (unsigned j = 1; j < vec.size(); j++) {
			delta.push_back(NU * point.d * (vec[j] - w[j]));
		}
		
		return delta;
	}

	void changeW(vector<vector<double>> delta) {
		vector<double> deltaRes = { 0.0, 0.0, 0.0 };

		for (vector<double> d : delta) {
			for (unsigned j = 0; j < d.size(); j++) {
				deltaRes[j] += d[j];
			}
		}

		for (unsigned j = 0; j < deltaRes.size(); j++) {
			w[j] += deltaRes[j] / delta.size();
		}

		normW();
	}

	int test(vector<double> vec) {
		vec = normalizeVector(vec);
		double cosPhi = 0.0;

		for (unsigned j = 1; j < vec.size(); j++) {
			cosPhi += w[j] * vec[j];
		}
		cout << cosPhi << endl;
		if (cosPhi > 0.95) {
			return 1;
		}

		return 0;
	}
};

int main() {
	srand((int)time(NULL));

	Grossberg neuron;

	vector<Point> ellipse1 = createEllipce(15.0, 7.5, 0.0, 10.0, 5.0, 1.0);
	vector<Point> ellipse2 = createEllipce(-5.0, 15.0, 0.0, 3.0, 12.0, 0.0);

	vector<vector<double>> delta;

	for (unsigned t = 0; t < N_ITER; t++) {
		for (Point p : ellipse1) {
			delta.push_back(neuron.teach(p));
		}
		for (Point p : ellipse2) {
			delta.push_back(neuron.teach(p));
		}
		
		neuron.changeW(delta);
	}
	
	vector<double> vec = neuron.getW();
	cout << "Vector W:" << endl;
	for (unsigned j = 1; j < vec.size(); j++) {
		cout << vec[j] << endl;
	}
	cout << endl;

	string answer;
	double x, y;
	bool questionFlag = true;

	while (true) {
		answer = "";

		if (questionFlag) {
			cout << "Do you want to test neuron? (y/n)" << endl;
			questionFlag = false;
		}
		getline(cin, answer);

		if (answer == "n") {
			break;
		}
		else if (answer == "y") {
			cout << "Enter x:" << endl;
			cin >> x;
			cout << "Enter y:" << endl;
			cin >> y;

			vector<double> v = { 0.0, x, y };

			if (neuron.test(v)) {
				cout << "Point belongs to the 1st group" << endl;
			}
			else {
				cout << "Point doesn't belong to the 1st group" << endl;
			}
			questionFlag = true;
		}
	}
}