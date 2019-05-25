//
// Created by lmunier on 06/05/19.
//

#ifndef PROJECT_VECTOR_HPP
#define PROJECT_VECTOR_HPP

#include <cstdio>
#include <string>
#include <iostream>

using namespace std;

template<typename Type>
class Vector2 {
public:
    Type x;
    Type y;

    explicit Vector2(Type x_val = 0.0f, Type y_val = 0.0f) :
            x(x_val),
            y(y_val) {}

    virtual void print() {
        printf("Vector2f : %f \t %f\n", this->x, this->y);
    }

#ifdef PRINT
    virtual string to_file() {
        string coord;
        coord = to_string(this->x) + "," + to_string(this->y) + "\n";

        return coord;
    }
#endif

    virtual Type norm() {
        return sqrtf(this->x*this->x + this->y*this->y);
    }

    Vector2<Type> operator-= (Vector2 vec) {
        this->x -= vec.x;
        this->y -= vec.y;

        return *this;
    }

    Vector2<Type> operator- (Vector2 vec) {
        Vector2 tmp = Vector2();

        tmp.x = this->x - vec.x;
        tmp.y = this->y - vec.y;

        return tmp;
    }

    Vector2<Type> operator+= (Vector2 vec) {
        this->x += vec.x;
        this->y += vec.y;

        return *this;
    }

    Vector2<Type> operator+ (Vector2 vec) {
        Vector2 tmp = Vector2();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;

        return tmp;
    }

    Vector2<Type> operator* (Type scalar) {
        Vector2 tmp = Vector2();

        tmp.x = scalar*this->x;
        tmp.y = scalar*this->y;

        return tmp;
    }
};

typedef Vector2<float> Vector2f;


template<typename Type>
class Vector3 : public Vector2<Type> {
public:
    Type z;

    explicit Vector3(Type x_val = 0.0f, Type y_val = 0.0f, Type z_val = 0.0f) :
        Vector2f(x_val, y_val),
        z(z_val) {}

    void print() override {
        printf("Vector3f : %f \t %f \t %f\n", this->x, this->y, this->z);
    }

#ifdef PRINT
    virtual string to_file() {
        string coord;
        coord = to_string(this->x) + "," + to_string(this->y) + "," + to_string(this->z) + "\n";

        return coord;
    }
#endif

    Type norm() override {
        return sqrtf(this->x*this->x + this->y*this->y + this->z*this->z);
    }

    Vector3<Type> operator-= (Vector3 vec) {
        this->x -= vec.x;
        this->y -= vec.y;
        this->z -= vec.z;

        return *this;
    }

    Vector3<Type> operator- (Vector3 vec) {
        Vector3 tmp = Vector3();

        tmp.x = this->x - vec.x;
        tmp.y = this->y - vec.y;
        tmp.z = this->z - vec.z;

        return tmp;
    }

    Vector3<Type> operator+= (Vector3 vec) {
        this->x += vec.x;
        this->y += vec.y;
        this->z += vec.z;

        return *this;
    }

    Vector3<Type> operator+ (Vector3 vec) {
        Vector3 tmp = Vector3();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;
        tmp.z = this->z + vec.z;

        return tmp;
    }

    Vector3<Type> operator* (Type scalar) {
        Vector3 tmp = Vector3();

        tmp.x = scalar*this->x;
        tmp.y = scalar*this->y;
        tmp.z = scalar*this->z;

        return tmp;
    }
};

typedef Vector3<float> Vector3f;

#endif //PROJECT_VECTOR_HPP
