#include <iostream>
#include <errno.h>
#include<iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stack>
#include <queue>
#include <cmath>
#define pi (acos(0.0)*2)
using namespace std;
// multi purpose class for point, vector in non homogenious way
struct Vector{
    string command;
    double x;
    double y;
    double z;
    double w;
};
struct Matrix{
    double cells[4][4];
};
struct Triangle{
    //just consider the x,y,z
    struct Vector point1;
    struct Vector point2;
    struct Vector point3;
    struct Matrix mat;
};

string s1 = "stage1.txt";
string s2= "stage2.txt";
string s3 = "stage3.txt";
ofstream file1(s1.c_str());
ofstream file2(s2.c_str());
ofstream file3(s3.c_str());
//gluLookAt params
struct Vector eye,look,up;
struct Vector l,r,u;
//gluPerspective params
double fovX,fovY,aspectRatio,near,far,t,rh;
//how many objects do I intend to draw
queue <struct Triangle> objects;
stack <int> index;
//transformation stack
stack <struct Matrix> operations;
//Identity Matrix
double I[4][4] = {{1,0,0,0},
                  {0,1,0,0},
                  {0,0,1,0},
                  {0,0,0,1}};
double V[4][4] = {{0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,1}};
double R[4][4] = {{0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,1}};
double T[4][4] = {{1,0,0,0},
                  {0,1,0,0},
                  {0,0,1,0},
                  {0,0,0,1}};
double P[4][4] = {{0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0}};
void printCoordinates(struct Vector point)
{
    cout<<point.x<<' '<<point.y<<' '<<point.z<<endl;
}

void printMatrix(double m[4][4])
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            cout<<m[i][j]<<' ';
        }
        cout<<endl;
    }
}
void crossProduct(struct Vector v3,struct Vector v1,struct Vector v2)
{
    v3.x = v1.y * v2.z - v1.z * v2.y;
    v3.y = v1.z * v2.x - v1.x * v2.z;
    v3.z = v1.x * v2.y - v1.y * v2.x;
}
double dotProduct(struct Vector v1,struct Vector v2)
{
    double res;
    res = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    return res;
}
void multVector(struct Vector v2,struct Vector v1,double m)
{
    v2.x = v1.x * m;
    v2.y = v1.y * m;
    v2.z = v1.z * m;
}
void subVector(struct Vector v3,struct Vector v1,struct Vector v2)
{
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    v3.z = v1.z - v2.z;
}

void pointToArray(Vector v,double point[4][1])
{
    point[0][0] = v.x;
    point[1][0] = v.y;
    point[2][0] = v.z;
    point[3][0] = v.w;
    //cout<<v.w<<endl;
}

void multiplyMatrix(double mult[4][4],double source[4][4],double dest[4][4])
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            mult[i][j]=0;
        }
    }
    //printMatrix(source);
    //printMatrix(dest);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                mult[i][j] += source[i][k] * dest[k][j];
            }
        }
    }
}
void multiplyTriangle(double mult[4][1],double dest[4][1],double source[4][4])
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 1; j++)
        {
            mult[i][j]=0;
        }
    }
    //printMatrix(source);
    //cout<<dest[0][0]<<' '<<dest[1][0]<<' '<<dest[2][0]<<' '<<dest[3][0]<<endl;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 1; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                mult[i][j] += source[i][k] * dest[k][j];
                //cout<<mult[i][j]<<' ';
            }
            //cout<<endl;
        }
    }
    //cout<<mult[0][0]<<' '<<mult[1][0]<<' '<<mult[2][0]<<' '<<mult[3][0]<<endl;
    //cout<<"print multiplied point"<<endl;
    //cout<<dest[0][0]<<' '<<dest[1][0]<<' '<<dest[3][0]<<endl;
}

void normalize(double point1[4][1],double point2[4][1])
{
    point2[0][0] = point1[0][0]/point1[3][0];
    point2[1][0] = point1[1][0]/point1[3][0];
    point2[2][0] = point1[2][0]/point1[3][0];
    point2[3][0] = point1[3][0]/point1[3][0];
}
void assignValues(double dest[4][4],double source[4][4])
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            dest[i][j] = source[i][j];
        }
    }
}
struct Vector Rot(struct Vector ax,struct Vector a,double angle)
{
    struct Vector res;
    float theta,cos_theta,inv,sin_theta;
    struct Vector cross;
    double dot,val;
    theta = (angle/180)*pi;
    cos_theta = cos(theta);
    if(abs(cos_theta) < 0.01)
        cos_theta = 0;
    inv = 1 - cos_theta;
    sin_theta = sin(theta);
    if(abs(sin_theta) < 0.01)
        sin_theta = 0;
    //cout<<"Angles"<<endl;
    //cout<<theta<<' '<<cos_theta<<' '<<sin_theta<<' '<<inv<<endl;
    //cout<<"Cross"<<endl;
    cross.x = a.y*ax.z - a.z*ax.y;
    cross.y = a.z*ax.x - a.x*ax.z;
    cross.z = a.x*ax.y - a.y*ax.x;
    //printCoordinates(cross);
    //val = sqrt(cross.x*cross.x+cross.y*cross.y+cross.z*cross.z);
    dot = a.x*ax.x + a.y*ax.y + a.z*ax.z;
    res.x = cos_theta*ax.x + inv*dot*a.x + sin_theta*cross.x;
    res.y = cos_theta*ax.y + inv*dot*a.y + sin_theta*cross.y;
    res.z = cos_theta*ax.z + inv*dot*a.z + sin_theta*cross.z;
    //val = sqrt(res.x*res.x+res.y*res.y+res.z*res.z);
    //res.x /= val;
    //res.y /= val;
    //res.z /= val;
    //cout<<"cross"<<endl;
    //printCoordinates(cross);
    return res;
}
void constructMatrix(double operation[4][4],struct Vector v)
{
    string op;
    op = v.command;
    if(op == "translate"){
        assignValues(operation,I);
        operation[0][3] = v.x;
        operation[1][3] = v.y;
        operation[2][3] = v.z;
        //operation[3][3] = v.w;
        //printMatrix(operation);
    }
    else if(op == "scale"){
        assignValues(operation,I);
        operation[0][0] = v.x;
        operation[1][1] = v.y;
        operation[2][2] = v.z;
        //operation[3][3] = v.w;
        //printMatrix(operation);
    }
    else if(op == "rotate"){
        assignValues(operation,I);
        struct Vector ax,a,c1,c2,c3;
        int val;
        double d;
        double theta = v.w;
        a.x = v.x/sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
        a.y = v.y/sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
        a.z = v.z/sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
        //cout<<a.x<<' '<<a.y<<' '<<a.z<<endl;
        //a.x = v.x;
        //a.y = v.y;
        //a.z = v.z;
        ax.x = 1;
        ax.y = 0;
        ax.z = 0;
        c1=Rot(ax,a,theta);
        operation[0][0] = c1.x;
        operation[1][0] = c1.y;
        operation[2][0] = c1.z;
        ax.x = 0;
        ax.y = 1;
        ax.z = 0;
        c2=Rot(ax,a,theta);
        operation[0][1] = c2.x;
        operation[1][1] = c2.y;
        operation[2][1] = c2.z;
        ax.x = 0;
        ax.y = 0;
        ax.z = 1;
        c3=Rot(ax,a,theta);
        operation[0][2] = c3.x;
        operation[1][2] = c3.y;
        operation[2][2] = c3.z;
        printMatrix(operation);
        //a_x = v.x;
        //a_y = v.y;
        //a_z = v.z;
        //cout<<pi<<endl;
    }
    else{
        cout<<"unrecognized op";
    }
}

void printObjects(struct Triangle t, struct Matrix m)
{
    //int itr = objects.size();
    //while(itr >0){
        //cout<<"Iteration "<<itr<<endl;
        //struct Triangle t;
        //t = objects.front();
        //struct Matrix m;
        //assignValues(m.cells,t.mat.cells);
        double p[4][1];
        double p11[4][1],p12[4][1],p13[4][1];
        double p21[4][1],p22[4][1],p23[4][1];
        double p31[4][1],p32[4][1],p33[4][1];
        cout<<"print stage 1 triangle"<<endl;
        pointToArray(t.point1,p);
        multiplyTriangle(p11,p,m.cells);
        cout<<p11[0][0]<<' '<<p11[1][0]<<' '<<p11[2][0]<<endl;
        file1<<p11[0][0]<<' '<<p11[1][0]<<' '<<p11[2][0]<<endl;
        pointToArray(t.point2,p);
        multiplyTriangle(p12,p,m.cells);
        cout<<p12[0][0]<<' '<<p12[1][0]<<' '<<p12[2][0]<<endl;
        file1<<p12[0][0]<<' '<<p12[1][0]<<' '<<p12[2][0]<<endl;
        pointToArray(t.point3,p);
        multiplyTriangle(p13,p,m.cells);
        cout<<p13[0][0]<<' '<<p13[1][0]<<' '<<p13[2][0]<<endl;
        file1<<p13[0][0]<<' '<<p13[1][0]<<' '<<p13[2][0]<<endl;
        cout<<endl;
        file1<<endl;
        //printMatrix(V);
        cout<<"print stage 2 triangle"<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p11,V);
        normalize(p,p21);
        cout<<p21[0][0]<<' '<<p21[1][0]<<' '<<p21[2][0]<<endl;
        file2<<p21[0][0]<<' '<<p21[1][0]<<' '<<p21[2][0]<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p12,V);
        normalize(p,p22);
        cout<<p22[0][0]<<' '<<p22[1][0]<<' '<<p22[2][0]<<endl;
        file2<<p22[0][0]<<' '<<p22[1][0]<<' '<<p22[2][0]<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p13,V);
        normalize(p,p23);
        cout<<p23[0][0]<<' '<<p23[1][0]<<' '<<p23[2][0]<<endl;
        file2<<p23[0][0]<<' '<<p23[1][0]<<' '<<p23[2][0]<<endl;
        cout<<endl;
        file2<<endl;
        cout<<"print stage 3 triangle"<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p21,P);
        normalize(p,p31);
        cout<<p31[0][0]<<' '<<p31[1][0]<<' '<<p31[2][0]<<endl;
        file3<<p31[0][0]<<' '<<p31[1][0]<<' '<<p31[2][0]<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p22,P);
        normalize(p,p32);
        cout<<p32[0][0]<<' '<<p32[1][0]<<' '<<p32[2][0]<<endl;
        file3<<p32[0][0]<<' '<<p32[1][0]<<' '<<p32[2][0]<<endl;
        p[0][0] = 0;
        p[1][0] = 0;
        p[2][0] = 0;
        p[3][0] = 0;
        multiplyTriangle(p,p23,P);
        normalize(p,p33);
        cout<<p33[0][0]<<' '<<p33[1][0]<<' '<<p33[2][0]<<endl;
        file3<<p33[0][0]<<' '<<p33[1][0]<<' '<<p33[2][0]<<endl;
        cout<<endl;
        file3<<endl;
        //objects.pop();
        //itr --;
    //}
}
void fileRead(string path)
{
    ifstream infile(path.c_str());
    string line;
    int line_count = 0;
    double a,b,c;
    string instruction;
    //reading first three lines for gluLookAt()
    while (line_count < 3)
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        switch(line_count)
        {
            case 0:
                eye.command = "gluLookAtParam-Eye";
                //eye.z = 0;
                iss >> eye.x >> eye.y >> eye.z;
                //printCoordinates(eye);
                T[0][3] = (-1) * eye.x;
                T[1][3] = (-1) * eye.y;
                T[2][3] = (-1) * eye.z;
                //printMatrix(T);
                break;
            case 1:
                look.command = "gluLookAtParam-Look";
                //look.z = 0;
                iss >> look.x >> look.y >> look.z;
                //printCoordinates(look);
                //subVector(l,look,eye);
                l.x = look.x - eye.x;
                l.y = look.y - eye.y;
                l.z = look.z - eye.z;
                //printCoordinates(l);
                //cout<<(sqrt((l.x*l.x)+(l.y*l.y)+(l.z*l.z)));
                double d;
                d = sqrt((l.x*l.x)+(l.y*l.y)+(l.z*l.z));
                l.x = l.x/d;
                l.y = l.y/d;
                    l.z = l.z/d;
                //printCoordinates(l);
                break;
            case 2:
                up.command = "gluLookAtParam-Up";
                //up.z = 0;
                iss >> up.x >> up.y >> up.z;
                //printCoordinates(up);
                //crossProduct(r,l,up);
                r.x = l.y*up.z - l.z*up.y;
                r.y = l.z*up.x - l.x*up.z;
                r.z = l.x*up.y - l.y*up.x;
                double e;
                e = sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
                r.x = r.x/e;
                r.y = r.y/e;
                r.z = r.z/e;
                //printCoordinates(r);
                //crossProduct(u,r,l);
                u.x = r.y*l.z - r.z*l.y;
                u.y = r.z*l.x - r.x*l.z;
                u.z = r.x*l.y - r.y*l.x;
                //printCoordinates(u);
                R[0][0] = r.x;
                R[0][1] = r.y;
                R[0][2] = r.z;
                R[0][3] = 0;
                R[1][0] = u.x;
                R[1][1] = u.y;
                R[1][2] = u.z;
                R[1][3] = 0;
                R[2][0] = -l.x;
                R[2][1] = -l.y;
                R[2][2] = -l.z;
                R[2][3] = 0;
                //printMatrix(T);
                //printMatrix(R);
                multiplyMatrix(V,R,T); //problem here
                //printMatrix(V);
                break;
            default:
                break;
        }
        //cout<<a<<' '<<b<<' '<<c<<endl;
        line_count ++;
    }
    //now we get gluPerspective() params
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> fovY >> aspectRatio >> near >>far;
    fovX = fovY * aspectRatio;
    t = near * tan(((double)fovY/(double)360)*pi);
    rh = near * tan(((double)fovX/(double)360)*pi);
    //cout<<t<<' '<<rh<<endl;
    P[0][0] = near/rh;
    P[1][1] = near/t;
    P[2][2] = -(far+near)/(far-near);
    P[2][3] = -(2*far*near)/(far-near);
    P[3][2] = -1;
    P[3][3] = 0;
    //printMatrix(P);
    //cout<<fovY<<' '<<aspectRatio<<' '<<near<<' '<<far<<endl;
    //Initialize Identity Matrix to the operations stack
    struct Matrix m;
    assignValues(m.cells,I);
    //printMatrix(m.cells);
    operations.push(m);
    //now we read All the 7 type of commands
    while(std::getline(infile, line))
    {
        std::istringstream iss(line);
        if(line == "triangle")
        {
            line_count = 0;
            struct Triangle t;
            while(line_count < 3)
            {
                std::getline(infile, line);
                std::istringstream iss(line);
                switch(line_count)
                {
                    case 0:
                        iss>> t.point1.x >> t.point1.y >> t.point1.z;
                        t.point1.w = 1;
                        //printCoordinates(t.point1);
                        break;
                    case 1:
                        iss>> t.point2.x >> t.point2.y >> t.point2.z;
                        t.point2.w = 1;
                        //printCoordinates(t.point2);
                        break;
                    case 2:
                        iss>> t.point3.x >> t.point3.y >> t.point3.z;
                        t.point3.w = 1;
                        //printCoordinates(t.point3);
                        break;
                    default:
                        break;
                }
                line_count ++;
            }
            struct Matrix m;
            m = operations.top();
            //assignValues(t.mat.cells,m.cells);
            //printMatrix(m.cells);
            printObjects(t,m);
            //printMatrix(t.mat.cells);
            objects.push(t);
        }
        else if(line == "translate")
        {
            std::getline(infile, line);
            std::istringstream iss(line);
            struct Vector v;
            struct Matrix m,n,p;
            n = operations.top();
            v.command = "translate";
            v.w = 1;
            iss >> v.x >> v.y >> v.z;
            constructMatrix(m.cells,v);
            //cout<<"Currently "<<endl;
            //printMatrix(n.cells);
            //cout<<"Just Translate "<<endl;
            //printMatrix(m.cells);
            multiplyMatrix(p.cells,n.cells,m.cells);
            //cout<<"After Multiplication "<<endl;
            //printMatrix(p.cells);
            operations.push(p);
        }
        else if(line == "scale")
        {
            std::getline(infile, line);
            std::istringstream iss(line);
            struct Vector v;
            struct Matrix m,n,p;
            n = operations.top();
            v.command = "scale";
            v.w = 1;
            iss >> v.x >> v.y >> v.z;
            constructMatrix(m.cells,v);
            //cout<<"Currently "<<endl;
            //printMatrix(n.cells);
            //cout<<"Just Scale "<<endl;
            //printMatrix(m.cells);
            multiplyMatrix(p.cells,n.cells,m.cells);
            //cout<<"After Multiplication "<<endl;
            //printMatrix(p.cells);
            operations.push(p);
        }
        else if(line == "rotate")
        {
            std::getline(infile, line);
            std::istringstream iss(line);
            struct Vector v;
            struct Matrix m,n,p;
            n = operations.top();
            v.command = "rotate";
            iss >> v.w >> v.x >> v.y >> v.z;
            constructMatrix(m.cells,v);
            //cout<<"Currently "<<endl;
            //printMatrix(n.cells);
            //cout<<"Just Rotate "<<endl;
            //printMatrix(m.cells);
            multiplyMatrix(p.cells,n.cells,m.cells);
            //cout<<"After Multiplication "<<endl;
            //printMatrix(p.cells);
            operations.push(p);
        }
        else if(line == "push")
        {
            int t;
            t = operations.size();
            //cout<<"push "<<t<<" current situation "<<endl;
            struct Matrix m = operations.top();
            //printMatrix(m.cells);
            index.push(t);
        }
        else if(line == "pop")
        {
            int t;
            t = index.top();
            //cout<<"pop "<<t<<endl;
            index.pop();
            while(operations.size() > t){
                operations.pop();
            }
            struct Matrix m;
            m = operations.top();
            //cout<<"After Pop"<<endl;
            //printMatrix(m.cells);
        }
        else if(line == "end")
        {

        }
        else
        {
            cout<<"unrecognized command"<<endl;
        }
    }
    infile.close();

}

int main()
{
    string path,opath1,opath2,opath3;
    path = "scene.txt";
    fileRead(path);
    file1.close();
    file2.close();
    file3.close();
    return 0;
}
