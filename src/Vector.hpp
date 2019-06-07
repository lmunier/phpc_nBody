/**
 * @file Vector.hpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
 *
 * 2D/3D Vector class to have a minimal implementation of vector with some functions.
 */

#ifndef PROJECT_VECTOR_HPP
#define PROJECT_VECTOR_HPP

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 */
#include "constants.hpp"

/**
 * @namespace std to simplify code implementation
 */
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
    Type v[DIM_2];

    /**
     * Surcharge constructor of the 2D vector. Call Vector2 constructor to construct it.
     *
     * @param val coordinate on all axis
     */
    explicit Vector2(Type val = 0.0f) {
        v[0] = val;
        v[1] = val;
    }

    /**
     * Surcharge  constructor of the 2D vector.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     */
    Vector2(Type x_val, Type y_val) {
        v[0] = x_val;
        v[1] = y_val;
    }

#ifdef PRINT
    /**
     * Print the coordinates of the given vector in the console.
     */
    void print() {
        printf("Vector2f : %f \t %f\n", this->v[0], this->v[1]);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    string to_file() {
        string coord;
        coord = to_string(this->v[0]) + "," + to_string(this->v[1]) + "\n";

        return coord;
    }
#endif

    /**
     * Compute norm of the given vector.
     *
     * @return norm of the given vector
     */
    Type norm() {
        return sqrtf(this->v[0]*this->v[0] + this->v[1]*this->v[1]);
    }

    Vector2<Type> pow(Type power) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = powf(this->v[0], power);
        tmp.v[1] = powf(this->v[1], power);

        return tmp;
    }

    /**
     * Override - to compute on vector2<Type>. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    Vector2<Type> operator- (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] - vec.v[0];
        if (abs(tmp.v[0]) < EPSILON)
            tmp.v[0] = 0.0f;

        tmp.v[1] = this->v[1] - vec.v[1];
        if (abs(tmp.v[1]) < EPSILON)
            tmp.v[1] = 0.0f;

        return tmp;
    }

    /**
     * Override + to compute on vector2<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    Vector2<Type> operator+ (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] + vec.v[0];
        tmp.v[1] = this->v[1] + vec.v[1];

        return tmp;
    }

    /**
     * Surcharge override + to compute on vector3<Type>.
     *
     * @param scalar Type to add to the given vector
     * @return given vector plus scalar argument on each dimension
     */
    Vector2<Type> operator+ (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] + scalar;
        tmp.v[1] = this->v[1] + scalar;

        return tmp;
    }

    /**
     * Override * to compute on vector2<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    Vector2<Type> operator* (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] * scalar;
        tmp.v[1] = this->v[1] * scalar;

        return tmp;
    }

    /**
     * Override / to compute on vector2<Type> with a scalar argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by scalar argument
     */
    Vector2<Type> operator/ (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] / scalar;
        tmp.v[1] = this->v[1] / scalar;

        return tmp;
    }

    /**
     * Override / to compute on vector2<Type> with a vector argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by vec argument
     */
    Vector2<Type> operator/ (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.v[0] = this->v[0] / vec.v[0];
        tmp.v[1] = this->v[1] / vec.v[1];

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
class Vector3 {
public:
    /**
     * Last coordinate of the vector in 3D.
     */
    Type v[DIM_3];

    /**
     * Surcharge constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param val coordinate on all axis
     */
    explicit Vector3(Type val = 0.0f) {
        v[0] = val;
        v[1] = val;
        v[2] = val;
    }

    /**
     * Surcharge constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     * @param y_val coordinate on z axis
     */
    Vector3(Type x_val, Type y_val, Type z_val) {
        v[0] = x_val;
        v[1] = y_val;
        v[2] = z_val;
    }

#ifdef PRINT
    /**
     * Print function to print value in console.
     */
    void print() {
        printf("Vector3f : %f \t %f \t %f\n", this->v[0], this->v[1], this->v[2]);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    string to_file() {
        string coord;
        coord = to_string(this->v[0]) + "," + to_string(this->v[1]) + "," + to_string(this->v[2]) + "\n";

        return coord;
    }
#endif
    /**
     * Norm computation of the given vector.
     *
     * @return norm of the given vector
     */
    Type norm() {
        return sqrtf(this->v[0]*this->v[0] + this->v[1]*this->v[1] + this->v[2]*this->v[2]);
    }

    Vector3<Type> pow(Type power) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = powf(this->v[0], power);
        tmp.v[1] = powf(this->v[1], power);
        tmp.v[2] = powf(this->v[2], power);

        return tmp;
    }

    /**
     * Override - to compute on vector3<Type    >. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    Vector3<Type> operator- (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = this->v[0] - vec.v[0];
        if (abs(tmp.v[0]) < EPSILON)
            tmp.v[0] = 0.0f;

        tmp.v[1] = this->v[1] - vec.v[1];
        if (abs(tmp.v[1]) < EPSILON)
            tmp.v[1] = 0.0f;

        tmp.v[2] = this->v[2] - vec.v[2];
        if (abs(tmp.v[2]) < EPSILON)
            tmp.v[2] = 0.0f;

        return tmp;
    }

    /**
     * Override + to compute on vector3<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    Vector3<Type> operator+ (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = this->v[0] + vec.v[0];
        tmp.v[1] = this->v[1] + vec.v[1];
        tmp.v[2] = this->v[2] + vec.v[2];

        return tmp;
    }

    /**
     * Surcharge override + to compute on vector3<Type>.
     *
     * @param scalar Type to add to the given vector
     * @return given vector plus scalar argument on each dimension
     */
    Vector3<Type> operator+ (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = this->v[0] + scalar;
        tmp.v[1] = this->v[1] + scalar;
        tmp.v[2] = this->v[2] + scalar;

        return tmp;
    }

    /**
     * Override * to compute on vector3<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    Vector3<Type> operator* (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = scalar * this->v[0];
        tmp.v[1] = scalar * this->v[1];
        tmp.v[2] = scalar * this->v[2];

        return tmp;
    }

    /**
     * Override / to compute on vector3<Type>.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by scalar argument
     */
    Vector3<Type> operator/ (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = this->v[0] / scalar;
        tmp.v[1] = this->v[1] / scalar;
        tmp.v[2] = this->v[2] / scalar;

        return tmp;
    }

    /**
     * Override / to compute on vector3<Type> with a vector argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by vec argument
     */
    Vector3<Type> operator/ (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.v[0] = this->v[0] / vec.v[0];
        tmp.v[1] = this->v[1] / vec.v[1];
        tmp.v[2] = this->v[2] / vec.v[2];

        return tmp;
    }
};

/**
 * @typedef define a typedef for 3D float vector
 */
typedef Vector3<float> Vector3f;

template<typename Type>
__global__ void add(Vector3<Type> *vec_i, Vector3<Type> *vec_j) {
    unsigned int i = threadIdx.x;
    printf("Coucou\n");
    printf("%f\n", vec_i->v[i]);

    printf("Hello\n");
    vec_i->v[i] += vec_j->v[i];
}

#endif //PROJECT_VECTOR_HPP
