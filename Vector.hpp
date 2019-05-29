/**
 * @file Vector.hpp
 * @author Munier Louis
 * @date 27.05.19
 * @version 1.0
 *
 * 2D/3D Vector class to have a minimal implementation of vector with some functions.
 */

#ifndef PROJECT_VECTOR_HPP
#define PROJECT_VECTOR_HPP

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 *
 * @namespace use std to simplify code implementation
 */
#include "constants.hpp"

using namespace std;

/**
 * Class of 2D vectors.
 *
 * @class Vector2
 * @tparam Type of the vector (int, float, etc ...)
 */
template<typename Type>
class Vector2 {
public:
    /**
     * coordinates of the vector.
     */
    Type x;
    Type y;

    /**
     * Constructor of the 2D vector.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     */
    explicit Vector2(Type x_val = 0.0f, Type y_val = 0.0f) :
            x(x_val),
            y(y_val) {}

#ifdef PRINT
    /**
     * Print the coordinates of the given vector in the console.
     */
    virtual void print() {
        printf("Vector2f : %f \t %f\n", this->x, this->y);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    virtual string to_file() {
        string coord;
        coord = to_string(this->x) + "," + to_string(this->y) + "\n";

        return coord;
    }
#endif

    /**
     * Compute norm of the given vector.
     *
     * @return norm of the given vector
     */
    virtual Type norm() {
        return sqrtf(this->x*this->x + this->y*this->y);
    }

    /**
     * Override - to compute on vector2<Type>. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    Vector2<Type> operator- (Vector2 vec) {
        Vector2 tmp = Vector2();

        tmp.x = this->x - vec.x;
        if (abs(tmp.x) < EPSILON)
            tmp.x = EPSILON;

        tmp.y = this->y - vec.y;
        if (abs(tmp.y) < EPSILON)
            tmp.y = EPSILON;

        return tmp;
    }

    /**
     * Override + to compute on vector2<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    Vector2<Type> operator+ (Vector2 vec) {
        Vector2 tmp = Vector2();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;

        return tmp;
    }

    /**
     * Override * to compute on vector2<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    Vector2<Type> operator* (Type scalar) {
        Vector2 tmp = Vector2();

        tmp.x = scalar*this->x;
        tmp.y = scalar*this->y;

        return tmp;
    }

    /**
     * Override / to compute on vector3<Type>.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by scalar argument
     */
    Vector2<Type> operator/ (Type scalar) {
        Vector2 tmp = Vector2();

        tmp.x = this->x/scalar;
        tmp.y = this->y/scalar;

        return tmp;
    }
};

/**
 * @typedef define a typedef for 2D float vector
 */
typedef Vector2<float> Vector2f;


/**
 * Class of 3D vectors.
 *
 * @class Vector3
 * @tparam Type of the vector (int, float, etc ...)
 */
template<typename Type>
class Vector3 : public Vector2<Type> {
public:
    /**
     * Last coordinate of the vector in 3D.
     */
    Type z;

    /**
     * Constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     * @param y_val coordinate on z axis
     */
    explicit Vector3(Type x_val = 0.0f, Type y_val = 0.0f, Type z_val = 0.0f) :
        Vector2f(x_val, y_val),
        z(z_val) {}

#ifdef PRINT
    /**
     * Override the print function of Vector2 to print 3D vectors.
     */
    void print() override {
        printf("Vector3f : %f \t %f \t %f\n", this->x, this->y, this->z);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    virtual string to_file() {
        string coord;
        coord = to_string(this->x) + "," + to_string(this->y) + "," + to_string(this->z) + "\n";

        return coord;
    }
#endif
/**
     * Override norm computation of the given vector.
     *
     * @return norm of the given vector
     */
    Type norm() override {
        return sqrtf(this->x*this->x + this->y*this->y + this->z*this->z);
    }

    /**
     * Override - to compute on vector3<Type>. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    Vector3<Type> operator- (Vector3 vec) {
        Vector3 tmp = Vector3();

        tmp.x = this->x - vec.x;
        if (abs(tmp.x) < EPSILON)
            tmp.x = EPSILON;

        tmp.y = this->y - vec.y;
        if (abs(tmp.y) < EPSILON)
            tmp.y = EPSILON;

        tmp.z = this->z - vec.z;
        if (abs(tmp.z) < EPSILON)
            tmp.z = EPSILON;

        return tmp;
    }

    /**
     * Override + to compute on vector3<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    Vector3<Type> operator+ (Vector3 vec) {
        Vector3 tmp = Vector3();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;
        tmp.z = this->z + vec.z;

        return tmp;
    }

    /**
     * Override * to compute on vector3<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    Vector3<Type> operator* (Type scalar) {
        Vector3 tmp = Vector3();

        tmp.x = scalar * this->x;
        tmp.y = scalar * this->y;
        tmp.z = scalar * this->z;

        return tmp;
    }

    /**
     * Override / to compute on vector3<Type>.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by scalar argument
     */
    Vector3<Type> operator/ (Type scalar) {
        Vector3 tmp = Vector3();

        tmp.x = this->x / scalar;
        tmp.y = this->y / scalar;
        tmp.z = this->z / scalar;

        return tmp;
    }
};

/**
 * @typedef define a typedef for 3D float vector
 */
typedef Vector3<float> Vector3f;

#endif //PROJECT_VECTOR_HPP
