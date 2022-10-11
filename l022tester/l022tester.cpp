// Riya Dev
// Period 5
// cd project2
// g++ -std=c++11 -o l022 -Wall l022.cpp
// ./l022

#include <iostream>
#include <fstream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>       /* pow */
#include <iomanip>      /* set precision*/
#include <climits>
#include <algorithm>    // std::max
#include <string>
#include <vector>
#include <list>

const int size = 800;
std::ofstream doc;

// point class
// private, getters, setters, constructors
class Point {
public:
    void setX(double x);
    double getX(void);
    void setY(double y);
    double getY(void);
    bool diff(Point p1, Point p2);
    Point(double x, double y);  // This is the constructor
    Point(void);
private:
    double x;
    double y;
};

Point::Point(void) {}

Point::Point(double xval, double yval) {
    x = xval;
    y = yval;
}

void Point::setX(double xval) {
    x = xval;
}
double Point::getX(void) {
    return x;
}
void Point::setY(double yval) {
    y = yval;
}
double Point::getY(void) {
    return y;
}

bool insidetriangle(double rx1, double ry1, double rx2, double ry2, double rx3, double ry3, double rx4, double ry4, double m23, double m34, double m24) {
    int count = 0;

    if (std::min(rx2, rx3) < rx1 && rx1 < std::max(rx2, rx3))
        if (rx1 < m23 * (rx1 - rx2) + ry2)
            count++;
    if (std::min(rx3, rx4) < rx1 && rx1 < std::max(rx3, rx4))
        if (m34 * (rx1 - rx3) + ry3 > rx1)
            count++;
    if (std::min(rx2, rx4) < rx1 && rx1 < std::max(rx2, rx4))
        if (m24 * (rx1 - rx4) + ry4 > rx1)
            count++;

    if (count == 3)
        return true;
    return false;
}

void part1() {
    std::ofstream file;
    file.open("points.txt");

    srand(time(NULL));
    double rx1 = double(rand()) / INT_MAX;
    double ry1 = double(rand()) / INT_MAX;

    double rx2 = double(rand()) / INT_MAX;
    double ry2 = double(rand()) / INT_MAX;

    double rx3 = double(rand()) / INT_MAX;
    double ry3 = double(rand()) / INT_MAX;

    double rx4 = double(rand()) / INT_MAX;
    double ry4 = double(rand()) / INT_MAX;

    double m12 = (ry1 - ry2) / (rx1 - rx2); //12
    double m23 = (ry2 - ry3) / (rx2 - rx3); //23
    double m13 = (ry1 - ry3) / (rx1 - rx3); //13
    double m14 = (ry1 - ry4) / (rx1 - rx4); //14
    double m24 = (ry2 - ry4) / (rx2 - rx4); //24
    double m34 = (ry3 - ry4) / (rx3 - rx4); //34

    //std::cout << rx1 << " " << ry1 << "  " << rx2 << " " << ry2 << "  " << rx3 << " " << ry3 << "  " << rx4 << " " << ry4 << "\n";

    //if a horizontal line intersects a slope one time, it is inside the triangle
    // triangle 234 - point 1
    std::cout << insidetriangle(rx1, ry1, rx2, ry2, rx3, ry3, rx4, ry4, m23, m34, m24) << "\n";
    while (insidetriangle(rx1, ry1, rx2, ry2, rx3, ry3, rx4, ry4, m23, m34, m24)) {
        rx1 = double(rand()) / INT_MAX;
        ry1 = double(rand()) / INT_MAX;
    }

    //triangle 134 - point 2
    std::cout << insidetriangle(rx2, ry2, rx1, ry1, rx3, ry3, ry4, ry4, m13, m34, m14) << "\n";
    while (insidetriangle(rx2, ry2, rx1, ry1, rx3, ry3, ry4, ry4, m13, m34, m14)) {
        rx2 = double(rand()) / INT_MAX;
        ry2 = double(rand()) / INT_MAX;
    }

    //triangle 124 - point 3
    std::cout << insidetriangle(rx3, ry3, rx1, ry1, rx2, ry2, rx4, ry4, m12, m24, m14) << "\n";
    while (insidetriangle(rx3, ry3, rx1, ry1, rx2, ry2, rx4, ry4, m12, m24, m14)) {
        rx3 = double(rand()) / INT_MAX;
        ry3 = double(rand()) / INT_MAX;
    }

    // triangle 123 - point 4
    std::cout << insidetriangle(rx4, ry4, rx1, ry1, rx2, ry2, rx3, ry3, m12, m23, m13) << "\n";
    while (insidetriangle(rx4, ry4, rx1, ry1, rx2, ry2, rx3, ry3, m12, m23, m13)) {
        rx4 = double(rand()) / INT_MAX;
        ry4 = double(rand()) / INT_MAX;
    }

    file << std::setprecision(17) << "(" << rx1 << "," << ry1 << ") , (" << rx2 << "," << ry2 << ") , (" << rx3 << "," << ry3 << ") , (" << rx4 << "," << ry4 << ")";
}

void illuminate(int** arr, int y, int x) {
    arr[x][y] = 0;
}

void notoutofbounds(int** array, int x, int y) {
    if (x >= 0 && x < size && y >= 0 && y < size)
        //array[x][y] = 0;
        illuminate(array, x, y);
}

// driving force x -> dx < dy

int** linecase1(int x1, int x2, int y1, int y2, int** arr) { // x1 < x2, y1 < y2, dx < dy = default bresenham's
    //std::cout << x1 << y1 << x2 << y2 << "\n";

    int dx = x2 - x1; // |
    int dy = y2 - y1; // _
    int e = dx - dy;
    int i = x1;

    for (int j = y1; j <= y2; j++) {
        notoutofbounds(arr, i, j);
        if (e >= 0) {
            i++;
            e -= dy;
        }
        e += dx;
    }
    return arr;
}

int** linecase2(int x1, int x2, int y1, int y2, int** arr) { // x1 < x2 && y1 > y2 && dx > dy
    //std::cout << x1 << y1 << x2 << y2 << "\n";

    int dx = x1 - x2; // |
    int dy = y2 - y1; // _
    int e = dx - dy;
    int i = x1;

    for (int j = y1; j <= y2; j++) {
        notoutofbounds(arr, i, j);
        if (e >= 0) {
            i--;
            e -= dy;
        }
        e += dx;
    }
    return arr;
}

int** linecase3(int x1, int x2, int y1, int y2, int** arr) { // y1 = y2 && x1 < x2 && dx < dy
    for (int i = x1; i < x2; i++) {
        //arr[i][y1] = 0;
        notoutofbounds(arr, i, y1);
    }
    return arr;
}

// driving force y -> dx > dy

int** linecase4(int x1, int x2, int y1, int y2, int** arr) { //x1 < x2 && y1 < y2 && dx > dy
    //std::cout << "case 4 " << x1 << y1 << " " << x2 << y2 << "\n";

    int dx = x2 - x1; // |
    int dy = y2 - y1; // _
    int e = dy - dx;
    int j = y1;

    for (int i = x1; i <= x2; i++) {
        //arr[i][j] = 0; // illuminate
        notoutofbounds(arr, i, j);
        if (e >= 0) {
            j++;
            e -= dx;
        }
        e += dy;
    }
    return arr;
}

int** linecase5(int x1, int x2, int y1, int y2, int** arr) {
    //std::cout << "case 5 " << x1 << y1 << " " << x2 << y2 << "\n";

    int dx = x2 - x1; // |
    int dy = y1 - y2; // _
    int e = dy - dx;
    int j = y1;

    for (int i = x1; i <= x2; i++) {
        //arr[i][j] = 0; // illuminate
        notoutofbounds(arr, i, j);
        if (e >= 0) {
            j--;
            e -= dx;
        }
        e += dy;
    }
    return arr;
}

int** linecase6(int x1, int x2, int y1, int y2, int** arr) { //x1 == x2 && y1 < y2
    //std::cout << "case 6" << x1 << y1 << " " << x2 << y2 << "\n";

    for (int j = y1; j < y2; j++) {
        //arr[x1][j] = 0;
        notoutofbounds(arr, x2, j);
    }
    return arr;
}

void chooselinefunction(int x1, int x2, int y1, int y2, int** array) {
    // the 6 freaking cases

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    // driving force x -> dx > dy

    // case 1
    if (x1 < x2 && y1 < y2 && dx <= dy) {
        //std::cout << "case 1a" << x1 << y1 << " " << x2 << y2 << "\n";
        array = linecase1(x1, x2, y1, y2, array);
    }
    else if (x2 < x1 && y2 < y1 && dx <= dy) {
        //std::cout << "case 1b" << x1 << y1 << " " << x2 << y2 << "\n";
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
        array = linecase1(x1, x2, y1, y2, array);
    }

    // case 2
    else if (x1 > x2 && y1 < y2 && dx <= dy) {
        //std::cout << "case 2a" << x1 << y1 << " " << x2 << y2 << "\n";
        array = linecase2(x1, x2, y1, y2, array);
    }
    else if (x1 < x2 && y1 > y2 && dx <= dy) {
        //std::cout << "case 2b" << x1 << y1 << " " << x2 << y2 << "\n";
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
        array = linecase2(x1, x2, y1, y2, array);
    }

    // case 3
    else if (y1 == y2 && x1 < x2) { // use higher precedence operator first
        //std::cout << "case 3a" << x1 << y1 << " " << x2 << y2 << "\n";
        array = linecase3(x1, x2, y1, y2, array);
    }
    else if (y1 == y2 && x1 > x2) {
        //std::cout << "case 3b" << x1 << y1 << " " << x2 << y2 << "\n";
        int temp = x1;
        x1 = x2;
        x2 = temp;
        array = linecase3(x1, x2, y1, y2, array);

    }

    // driving force y -> dy > dx

    // case 4
    else if (x1 < x2 && y1 < y2 && dx > dy) {
        //std::cout << "case 4a" << x1 << y1 << " " << x2 << y2 << "\n";
        array = linecase4(x1, x2, y1, y2, array);
    }
    else if (x2 < x1 && y2 < y1 && dx > dy) {
        //std::cout << "case 4b" << x1 << y1 << " " << x2 << y2 << "\n";
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
        array = linecase4(x1, x2, y1, y2, array);
    }

    // case 5
    else if (x1 < x2 && y1 > y2 && dx > dy) {
        //std::cout << "case 5a" << x1 << y1 << " " << x2 << y2 << "\n";
        array = linecase5(x1, x2, y1, y2, array);
    }
    else if (x2 < x1 && y2 > y1 && dx > dy) {
        //std::cout << "case 5b" << x1 << y1 << " " << x2 << y2 << "\n";
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
        array = linecase5(x1, x2, y1, y2, array);
    }

    // case 6
    else if (x1 == x2 && y1 < y2) { // use higher precedence operator first
        array = linecase6(x1, x2, y1, y2, array);
    }
    else if (x1 == x2 && y2 < y1) {
        int temp = y1;
        y1 = y2;
        y2 = temp;
        array = linecase6(x1, x2, y1, y2, array);
    }
}

void eulerline(int** array, int x1, int y1, int x2, int y2) {
    //centroid and circumcenter coordiinates
    //x1, x2, y1, y2, array

    double m = (double)(y1 - y2) / (x1 - x2);
    int yi = m * -x1 + y1;
    int yf = m * (size - x1 - 1) + y1;
    //std::cout << yi << " " << yf;
    // y - y1 = m (x - x1) -> y = mx - mx1 + y1
    // y - y2 = m(x - x2)  -> y = mx - mx2 + y2

    chooselinefunction(0, size - 1, yi, yf, array);
}

double distance(double x1, double y1, double x2, double y2) {
    double t1 = pow(x1 - x2, 2);
    double t2 = pow(y1 - y2, 2);
    //std::cout << t1 << " " << t2 << " " << sqrt(t1 + t2) << "\n";
    return (double)(sqrt(t1 + t2));
}

void circle(int** array, int xcenter, int ycenter, int R) {
    int y = R;
    int xmax = (int)(R * 0.70710678);
    int ysquare = y * y;
    int ty = (2 * y) - 1;
    int y2_new = ysquare;

    //std::cout << y << " " << xmax << " " << ysquare << " " << ty << " " << y2_new << "\n";

    for (int x = 0; x <= xmax + 1; x++) {
        if ((ysquare - y2_new) >= ty) {
            ysquare -= ty;
            y -= 1;
            ty -= 2;
        }
        notoutofbounds(array, y + xcenter, x + ycenter);
        notoutofbounds(array, -y + xcenter, x + ycenter);
        notoutofbounds(array, y + xcenter, -x + ycenter);
        notoutofbounds(array, -y + xcenter, -x + ycenter);

        notoutofbounds(array, x + xcenter, y + ycenter);
        notoutofbounds(array, -x + xcenter, y + ycenter);
        notoutofbounds(array, x + xcenter, -y + ycenter);
        notoutofbounds(array, -x + xcenter, -y + ycenter);

        y2_new -= (2 * x) - 3;
    }
}

std::vector<double> makesquarehelper(int** array, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double Ex, double Ey) {
    // step 4 
    // unite the unused point (D) with E and this is the 1st line
    double mDE = (double)(Dy - Ey) / (Dx - Ex);
    //eulerline(array, (int)(Dx * size), (int)(Dy * size), (int)(Ex * size), (int)(Ey * size));
        // y - Dy = mDE(x - Dx)
    //eulerline(array, (int)(Bx * size), (int)(By * size), (int)((1 + Bx) * size), (int)((mDE + By) * size));
        // y - By = mDE(x - Bx)

    // step 5
    // from A & C (step 2) draw perpendicular on DE
    double perpDE = -(1 / mDE);
    //eulerline(array, (int)(Ax * size), (int)(Ay * size), (int)((1 + Ax) * size), (int)((perpDE + Ay) * size));
          // y - Ay = perpDE(x - Ax)
    //eulerline(array, (int)(Cx * size), (int)(Cy * size), (int)((1 + Cx) * size), (int)((perpDE + Cy) * size));
        // y - Cy = perpDE(x - Cx)

    // step 6
    // find area using lines AB and AD
    double corner1x = (double)(Ay - perpDE * Ax - By + mDE * Bx) / (mDE - perpDE);  // AB
    double corner1y = By + mDE * (corner1x - Bx);
    Point corner1(corner1x, corner1y);
    //circle(array, (int)(corner1x * size), (int)(corner1y * size), 8);

    double corner2x = (double)(Ay - perpDE * Ax - Dy + mDE * Dx) / (mDE - perpDE);  // AD
    double corner2y = Dy + mDE * (corner2x - Dx);
    Point corner2(corner2x, corner2y);
    //circle(array, (int)(corner2x * size), (int)(corner2y * size), 8);
    
    double corner3x = (double)(Cy - perpDE * Cx - Dy + mDE * Dx) / (mDE - perpDE);  // CD
    double corner3y = Cy + perpDE * (corner3x - Cx);
    Point corner3(corner3x, corner3y);
    //circle(array, (int)(corner3x * size), (int)(corner3y * size), 8);
    
    double corner4x = (double)(Cy - perpDE * Cx - By + mDE * Bx) / (mDE - perpDE);  // CB
    double corner4y = Cy + perpDE * (corner4x - Cx);
    Point corner4(corner4x, corner4y);
    //circle(array, (int)(corner4x * size), (int)(corner4y * size), 8);

    double dist = distance(corner1x, corner1y, corner2x, corner2y);

    std::vector<double> vec;
    vec = {dist, corner1x, corner1y, corner2x, corner2y, corner3x, corner3y, corner4x, corner4y};
    doc << std::setprecision(17) << "(" << vec[1] << "," << vec[2] << ") , (" << vec[3] << "," << vec[4] << ") , (" << vec[5] << "," << vec[6] << ") , (" << vec[7] << "," << vec[8] << ") Area=" << pow(vec[0], 2) << "\n";
    return vec;
}

std::vector<double> makesquares(int** array, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy) {
    // step 1
    // pick 2 points (done by method header)
    double mAC = (double)(Cy - Ay)/(Cx - Ax);
    double lenAC = distance(Ax, Ay, Cx, Cy);

    // step 2
    // unite the two poitns (AC) and pick one of the other (B or D, doesn't matter) and draw a perpendicular line onto AC
    double mBE = -1 / mAC;
    
    // step 3
    // pick E on the line in step 2  such that |BE| = |AC|
    double dx = (double)lenAC/(sqrt(1 + pow(mBE, 2)));
    double dy = (double)(lenAC * mBE) / (sqrt(1 + pow(mBE, 2)));
    double Ex1;
    double Ex2;
    // up
    if (mBE < 0)
        Ex1 = Bx - dx;
    else
        Ex1 = Bx + dx;
    double Ey1 = By + dy;

    // down
    if (mBE < 0)
        Ex2 = Bx + dx;
    else
        Ex2 = Bx - dx;
    double Ey2 = By - dy;

    //std::cout << "dx " << dx << " dy " << dy << " Ex1 " << Ex1 << " Ey1 " << Ey1 << "\n";

    //circle(array, (int)(Ex1 * size), (int)(Ey1 * size), 4);
    //eulerline(array, (int)(Bx * size), (int)(By * size), (int)(Ex1 * size), (int)(Ey1 * size));
    //circle(array, (int)(Ex2 * size), (int)(Ey2 * size), 4);

    std::vector<double> topC = makesquarehelper(array, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex1, Ey1);
    std::vector<double> bottomC = makesquarehelper(array, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex2, Ey2);

    if (topC.front() < bottomC.front())
        return topC;
    return bottomC;
    //std::cout << topC << "\n";

    /*
    // diagonal 1
    // circle 1
    int xc12 = (int)(x1 * size + x2 * size) / 2;
    int yc12 = (int)(y1 * size + y2 * size) / 2;
    int r12 = distance(xc12, yc12, x1 * size, y1 * size);
    //circle(array, xc12, yc12, r12);

    double m12 = (double)(y2 - y1) / (x2 - x1);
    m12 = -1 / m12;
    int yf12 = ((r12 * m12) / sqrt(1 + pow(m12, 2))) + yc12;
    int xf12 = ((yf12 - yc12) / m12) + xc12;
    //circle(array, xf12, yf12, 4);

    // circle 2
    int xc34 = (int)(x3 * size + x4 * size) / 2;
    int yc34 = (int)(y3 * size + y4 * size) / 2;
    int r34 = distance(xc34, yc34, x3 * size, y3 * size);
    //circle(array, xc34, yc34, r34);

    double m34 = (double)(y4 - y3) / (x4 - x3);
    m34 = -1 / m34;
    int yf34 = ((r34 * m34) / sqrt(1 + pow(m34, 2))) + yc34;
    int xf34 = ((yf34 - yc34) / m34) + xc34;
    //circle(array, xf34, yf34, 4);

    // chooselinefunction(xf12, xf34, yf12, yf34, array);
    //eulerline(array, xf12, yf12, xf34, yf34);

    // diagonal 2
    int xc = (xf12 + xf34) / 2;
    int yc = (yf12 + yf34) / 2;

    double m = (double)(yf12 - yf34) / (xf12 - xf34);
    m = -1 / m;
    double d = distance(xf12, yf12, xf34, yf34);
    d = d / 2;
    // up
    int dy = (m * d) / (sqrt(1 + pow(m, 2)));
    int yf56 = dy + yc;
    int xf56 = ((yf56 - yc) / m) + xc;
    //circle(array, xf56, yf56, 4);
    int yf78 = yc - dy;
    int xf78 = ((yf78 - yc) / m) + xc;
    //circle(array, xf78, yf78, 4);
    //chooselinefunction(xf56, xf78, yf56, yf78, array);

    eulerline(array, xf12, yf12, xf56, yf56);
    eulerline(array, xf56, yf56, xf34, yf34);
    eulerline(array, xf34, yf34, xf78, yf78);
    eulerline(array, xf78, yf78, xf12, yf12);
    double imtired = distance(xf12, yf12, xf34, yf34);
    */
}

void wedrawing(int** array, double corner1x, double corner1y, double corner2x, double corner2y, double corner3x, double corner3y, double corner4x, double corner4y) {    
    eulerline(array, (int)(corner1x * size), (int)(corner1y * size), (int)(corner2x * size), (int)(corner2y * size));
    eulerline(array, (int)(corner2x * size), (int)(corner2y * size), (int)(corner3x * size), (int)(corner3y * size));
    eulerline(array, (int)(corner3x * size), (int)(corner3y * size), (int)(corner4x * size), (int)(corner4y * size));
    eulerline(array, (int)(corner4x * size), (int)(corner4y * size), (int)(corner1x * size), (int)(corner1y * size));
}

void part2() {
    //pulling the numbers from the file

    std::vector<std::string> vec;
    std::string line;
    std::ifstream points("points.txt");
    if (points.is_open()) {
        while (getline(points, line, ',')) {
            //std::cout << line << '\n';
            vec.push_back(line);
        }
        points.close();
    }
    else
        std::cout << "Unable to open file";

    doc.open("output.txt");
    doc << std::fixed;

    double x1 = std::stod(vec[0].substr(1, 19));
    double y1 = std::stod(vec[1].substr(0, 19));
    Point A(x1, y1);
    double x2 = std::stod(vec[2].substr(3, 19));
    double y2 = std::stod(vec[3].substr(0, 19));
    Point B(x2, y2);
    double x3 = std::stod(vec[4].substr(3, 19));
    double y3 = std::stod(vec[5].substr(0, 19));
    Point C(x3, y3);
    double x4 = std::stod(vec[6].substr(3, 19));
    double y4 = std::stod(vec[7].substr(0, 19));
    Point D(x4, y4);

    doc << std::setprecision(17) << "(" << x1 << "," << y1 << ") , (" << x2 << "," << y2 << ") , (" << x3 << "," << y3 << ") , (" << x4 << "," << y4 << ")\n";

    // ppm file
    std::ofstream file;
    file.open("output.ppm");
    file << "P3 " << size << " " << size << " 1\n";

    int** array;
    array = new int* [size];

    // fill with 1s
    for (int row = 0; row < size; ++row)
        array[row] = new int[size];
    for (int row = 0; row < size; ++row)
        for (int col = 0; col < size; ++col)
            array[row][col] = 1;

    // making circles for the 4 points
    circle(array, x1 * size, y1 * size, 2);
    circle(array, x2 * size, y2 * size, 2);
    circle(array, x3 * size, y3 * size, 2);
    circle(array, x4 * size, y4 * size, 2);

    //std::cout << "(" << x1 << ", " << y1 << ")" << "(" << x2 << ", " << y2 << ")" << "(" << x3 << ", " << y3 << ")" << "(" << x4 << ", " << y4 << ")" << "\n";

    // makesquares(array, x1, y1, x2, y2, x3, y3, x4, y4);
    std::vector<double> case1 = makesquares(array, A.getX(), A.getY(), B.getX(), B.getY(), C.getX(), C.getY(), D.getX(), D.getY());
    std::vector<double> case2 = makesquares(array, A.getX(), A.getY(), C.getX(), C.getY(), B.getX(), B.getY(), D.getX(), D.getY());
    std::vector<double> case3 = makesquares(array, A.getX(), A.getY(), B.getX(), B.getY(), D.getX(), D.getY(), C.getX(), C.getY());

    //std::vector<double> final;


    if (case1.front() < case2.front() && case1.front() < case3.front()) {
        wedrawing(array, case1[1], case1[2], case1[3], case1[4], case1[5], case1[6], case1[7], case1[8]);
        //final = case1;
    }
    else if (case2.front() < case1.front() && case2.front() < case3.front()) {
        wedrawing(array, case2[1], case2[2], case2[3], case2[4], case2[5], case2[6], case2[7], case2[8]);
        //final = case2;
    }
    else {
        wedrawing(array, case3[1], case3[2], case3[3], case3[4], case3[5], case3[6], case3[7], case3[8]);
        //final = case3;
    }

    // print array in output.ppm
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j)
            file << array[i][j] << " " << array[i][j] << " " << array[i][j] << ' ';
        file << std::endl;
    }
    file << "\n";
    file.close();

    doc.close();
}


int main() {
    //part1();
    part2();
}