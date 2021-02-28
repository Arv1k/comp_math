#include <iostream>
#include <cmath>
#include <vector>

// I don't like, but it's easier with it
using namespace std;

struct Pendulum {
private:
    // the parameters of the equation (l, L (lambda), g)
    // l*u''(t) + 2L*u'(t) + g*sin(u(t)) = 0
    double l = NAN;
    double L = NAN;
    double g = NAN;

    // the specified parameters (l -> 1, L/l -> y (gamma), g/l -> w (omega))
    // u'(t) = a(t)
    // a'(t) + 2y*a(t) + w*sin(u(t)) = 0
    double y;
    double w;

    // the initial conditions (u(0) = u0, u'(0) = a(0) = v0)
    double u0 = NAN;
    double v0 = NAN;

    // the period (T) and sampling step (ss)
    double T;
    double ss;

public:
    Pendulum();
    ~Pendulum();

    void ExplicitEuler(vector<double>& u) const;
    void RungeKutta(vector<double>& u) const;
    void AccurateSolution(vector<double>& u) const;

    // just a crutch
    double getT()  const { return T; };
    double getSS() const { return ss; };

    void showParams() const;
    void visualization(vector<double> dataset) const;
};

Pendulum::Pendulum() {
    cout << "Please enter the parameters ";
    cout << "(l, L (lambda), g, u0, v0): " << std::endl;
    cin >> l >> L >> g >> u0 >> v0;

    y  = L/l;
    w  = g/l;
    T  = 2*M_PI*sqrt(l/g);
    ss = 3*T/60;
}

Pendulum::~Pendulum() {
    l = NAN;
    L = NAN;
    g = NAN;

    u0 = NAN;
    v0 = NAN;

    y  = NAN;
    w  = NAN;

    T  = NAN;
    ss = NAN;
}

void Pendulum::ExplicitEuler(vector<double>& u) const {
    double un = u0;
    double an = v0;

    for (double t = 0; t < 3*T; t += ss) {
        u.push_back(un);

        an = an*(1 - ss*2*y) - ss*w*sin(un);
        un = un + ss*an;
    }
}

void Pendulum::RungeKutta(vector<double>& u) const {
    double un = u0;
    double an = v0;

    for (double t = 0; t < 3*T; t += ss) {
        u.push_back(un);

        double k1_a = -(2*y*an + w*sin(un));
        double k1_u = an;

        double k2_a = -(2*y*(an + ss/2*k1_a) + w*sin(un + ss/2*k1_u));
        double k2_u = an + ss/2*k1_a;

        double k3_a = -(2*y*(an + ss/2*k2_a) + w*sin(un + ss/2*k2_u));
        double k3_u = an + ss/2*k2_a;

        double k4_a = -(2*y*(an + ss*k3_a) + w*sin(un + ss*k3_u));
        double k4_u = an + ss*k3_a;

        an += ss/6*(k1_a + 2*k2_a + 2*k3_a + k4_a);
        un += ss/6*(k1_u + 2*k2_u + 2*k3_u + k4_u);
    }
}

/**
 * because of the doubles,
 * I use a handmade comparison method
 * (I consider the numbers equal if their difference
 * fits into a certain interval
 */
void Pendulum::AccurateSolution(vector<double>& u) const {
    if (y - sqrt(w) < -0.00005) {
        // W (big omega)
        double W = sqrt(w - pow(y, 2));

        double C1 = u0;
        double C2 = (v0 + u0*y)/W;

        for (double t = 0; t < 3*T; t += ss) {
            u.push_back((C1*cos(W*t) + C2*sin(W*t))*exp(-y*t));
        }
    } else if (y - sqrt(w) > 0.00005) {
        // Y (big gamma)
        double Y = sqrt(pow(y, 2) - w);

        double C1 = (-v0 + (Y - y)*u0)/(2*Y);
        double C2 = (w + u0*(y + Y))/(2*Y);

        for (double t = 0; t < 3*T; t += ss) {
           u.push_back((C1*exp(-Y*t) + C2*exp(Y*t))*exp(-y*t));
        }
    } else {
        double C1 = v0 + u0*y;
        double C2 = u0;

        for (double t = 0; t < 3*T; t += ss) {
            u.push_back((C1*t + C2)*exp(-y*t));
        }
    }
}

/**
 * just for self-check
 */
void Pendulum::showParams() const {
    cout << "Parameters:" << endl;
    cout << "l = " << l << endl;
    cout << "L = " << L << endl;
    cout << "g = " << g << endl << endl;

    cout << "u0 = " << u0 << endl;
    cout << "v0 = " << v0 << endl << endl;

    cout << "y = " << y << endl;
    cout << "w = " << w << endl;
    cout << "T = " << T << endl;
    cout << "ss = " << ss << endl << endl;
}

/**
 * I use my friend's function (for self-check too)
 */
void Pendulum::visualization(vector<double> dataset) const {
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

    fprintf(gnuplotPipe, "set multiplot layout 1, 1\n");
    fprintf(gnuplotPipe, "plot '-'\n");

    for (unsigned int i = 0; i < dataset.size(); i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", i * ss, dataset[i]);
    }

    fprintf(gnuplotPipe, "e\n");
    pclose(gnuplotPipe);
}

int main() {
    Pendulum ticker;
    ticker.showParams();

    vector<double> EE;
    ticker.ExplicitEuler(EE);

    vector<double> RK;
    ticker.RungeKutta(RK);

    vector<double> AS;
    ticker.AccurateSolution(AS);

    int precision = 4;
    cout.width(precision*3);
    cout << left << "t";
    cout.width(precision*3);
    cout<< left << "EE";
    cout.width(precision*3);
    cout<< left << "RK";
    cout.width(precision*3);
    cout<< left << "AS" << endl;

    for (size_t i = 0; i < AS.size(); i++) {
        cout.width(precision*3);
        cout.precision(precision);
        cout << left << ticker.getSS()*i;

        cout.width(precision*3);
        cout.precision(precision);
        cout << left << EE[i];

        cout.width(precision*3);
        cout.precision(precision);
        cout << left << RK[i];

        cout.width(precision*3);
        cout.precision(precision);
        cout << left << AS[i] << endl;
    }

    /*ticker.visualization(EE);
    ticker.visualization(RK);
    ticker.visualization(AS);*/
}
