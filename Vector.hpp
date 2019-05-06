//
// Created by lmunier on 06/05/19.
//

#ifndef PROJECT_VECTOR_HPP
#define PROJECT_VECTOR_HPP

#include <cstdio>


template<typename Type>
class Vector2 {
public:
    float x;
    float y;

    explicit Vector2(float x_val = 0.0f, float y_val = 0.0f) :
            x(x_val),
            y(y_val) {}

    virtual void print() {
        printf("Vector2f : %f \t %f", this->x, this->y);
    }

    Vector2<float> operator-= (Vector2 vec) {
        this->x -= vec.x;
        this->y -= vec.y;

        return *this;
    }

    Vector2<float> operator- (Vector2 vec) {
        Vector2 tmp = Vector2(0.0f, 0.0f);

        tmp.x = this->x - vec.x;
        tmp.y = this->y - vec.y;

        return tmp;
    }

private:
};

typedef Vector2<float> Vector2f;


template<typename Type>
class Vector3 : public Vector2<Type> {
public:
    float z;

    explicit Vector3(float x_val = 0.0f, float y_val = 0.0f, float z_val = 0.0f) :
        Vector2f(x_val, y_val),
        z(z_val) {}

    void print() override {
        printf("Vector3f : %f \t %f \t %f", this->x, this->y, this->z);
    }

    Vector3<float> operator-= (Vector3 vec) {
        this->x -= vec.x;
        this->y -= vec.y;
        this->z -= vec.z;

        return *this;
    }

    Vector3<float> operator- (Vector3 vec) {
        Vector3 tmp = Vector3(0.0f, 0.0f, 0.0f);

        tmp.x = this->x - vec.x;
        tmp.y = this->y - vec.y;
        tmp.z = this->z - vec.z;

        return tmp;
    }

private:
};

typedef Vector3<float> Vector3f;

#endif //PROJECT_VECTOR_HPP
