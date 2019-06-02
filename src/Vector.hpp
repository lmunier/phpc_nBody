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
    Type x;
    Type y;

    /**
     * Surcharge constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param val coordinate on all axis
     */
    __host__ __device__ explicit Vector2(Type val = 0.0f) : x(val), y(val) {}

    /**
     * Surcharge  constructor of the 2D vector.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     */
    __host__ __device__ Vector2(Type x_val, Type y_val) :
            x(x_val),
            y(y_val) {}

    __host__ void *operator new(size_t len) {
        void *ptr;
        cudaMallocManaged(&ptr, len);
        cudaDeviceSynchronize();
        return ptr;
    }

    __host__ void operator delete(void *ptr) {
        cudaDeviceSynchronize();
        cudaFree(ptr);
    }

#ifdef PRINT
    /**
     * Print the coordinates of the given vector in the console.
     */
    __host__ __device__ virtual void print() {
        printf("Vector2f : %f \t %f\n", this->x, this->y);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    __host__ __device__ virtual string to_file() {
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
    __device__ virtual Type norm() {
        return sqrtf(this->x*this->x + this->y*this->y);
    }

    __host__ __device__ Vector2<Type> pow(Type power) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = powf(this->x, power);
        tmp.y = powf(this->y, power);

        return tmp;
    }

    /**
     * Override - to compute on vector2<Type>. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    __host__ __device__ Vector2<Type> operator- (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = this->x - vec.x;
        if (abs(tmp.x) < EPSILON)
            tmp.x = 0.0f;

        tmp.y = this->y - vec.y;
        if (abs(tmp.y) < EPSILON)
            tmp.y = 0.0f;

        return tmp;
    }

    /**
     * Override + to compute on vector2<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    __host__ __device__ Vector2<Type> operator+ (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;

        return tmp;
    }

    /**
     * Surcharge override + to compute on vector3<Type>.
     *
     * @param scalar Type to add to the given vector
     * @return given vector plus scalar argument on each dimension
     */
    __host__ __device__ Vector2<Type> operator+ (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = this->x + scalar;
        tmp.y = this->y + scalar;

        return tmp;
    }

    /**
     * Override * to compute on vector2<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    __host__ __device__ Vector2<Type> operator* (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = scalar * this->x;
        tmp.y = scalar * this->y;

        return tmp;
    }

    /**
     * Override / to compute on vector2<Type> with a scalar argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by scalar argument
     */
    __host__ __device__ Vector2<Type> operator/ (Type scalar) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = this->x / scalar;
        tmp.y = this->y / scalar;

        return tmp;
    }

    /**
     * Override / to compute on vector2<Type> with a vector argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by vec argument
     */
    __host__ __device__ Vector2<Type> operator/ (Vector2<Type> vec) {
        Vector2<Type> tmp = Vector2<Type>();

        tmp.x = this->x / vec.x;
        tmp.y = this->y / vec.y;

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
     * Surcharge constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param val coordinate on all axis
     */
    __host__ __device__ explicit Vector3(Type val = 0.0f) :
            Vector2f(val, val),
            z(val) {}

    /**
     * Surcharge constructor of the 3D vector. Call Vector2 constructor to construct it.
     *
     * @param x_val coordinate on x axis
     * @param y_val coordinate on y axis
     * @param y_val coordinate on z axis
     */
    __host__ __device__ Vector3(Type x_val, Type y_val, Type z_val) :
        Vector2f(x_val, y_val),
        z(z_val) {}

#ifdef PRINT
    /**
     * Override the print function of Vector2 to print 3D vectors.
     */
    __host__ __device__ void print() override {
        printf("Vector3f : %f \t %f \t %f\n", this->x, this->y, this->z);
    }

    /**
     * Save coordinates of the given vector in a csv file.
     *
     * @return string of the given vector's coordinates
     */
    __host__ __device__ string to_file() override {
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
    __device__ Type norm() override {
        return sqrtf(this->x*this->x + this->y*this->y + this->z*this->z);
    }

    __device__ Vector3<Type> pow(Type power) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = powf(this->x, power);
        tmp.y = powf(this->y, power);
        tmp.z = powf(this->z, power);

        return tmp;
    }

    /**
     * Override - to compute on vector3<Type>. Avoid having a too small vector which reach to segfault.
     *
     * @param vec vector to subtract to the given vector
     * @return given vector minus vec argument
     */
    __host__ __device__ Vector3<Type> operator- (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = this->x - vec.x;
        if (abs(tmp.x) < EPSILON)
            tmp.x = 0.0f;

        tmp.y = this->y - vec.y;
        if (abs(tmp.y) < EPSILON)
            tmp.y = 0.0f;

        tmp.z = this->z - vec.z;
        if (abs(tmp.z) < EPSILON)
            tmp.z = 0.0f;

        return tmp;
    }

    /**
     * Override + to compute on vector3<Type>.
     *
     * @param vec vector to add to the given vector
     * @return given vector plus vec argument
     */
    __host__ __device__ Vector3<Type> operator+ (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = this->x + vec.x;
        tmp.y = this->y + vec.y;
        tmp.z = this->z + vec.z;

        return tmp;
    }

    /**
     * Surcharge override + to compute on vector3<Type>.
     *
     * @param scalar Type to add to the given vector
     * @return given vector plus scalar argument on each dimension
     */
    __host__ __device__ Vector3<Type> operator+ (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = this->x + scalar;
        tmp.y = this->y + scalar;
        tmp.z = this->z + scalar;

        return tmp;
    }

    /**
     * Override * to compute on vector3<Type>.
     *
     * @param Type value to multiply the given vector
     * @return given vector multiply by scalar argument
     */
    __host__ __device__ Vector3<Type> operator* (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

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
    __host__ __device__ Vector3<Type> operator/ (Type scalar) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = this->x / scalar;
        tmp.y = this->y / scalar;
        tmp.z = this->z / scalar;

        return tmp;
    }

    /**
     * Override / to compute on vector3<Type> with a vector argument.
     *
     * @param Type value to divide the given vector
     * @return given vector divide by vec argument
     */
    __host__ __device__ Vector3<Type> operator/ (Vector3<Type> vec) {
        Vector3<Type> tmp = Vector3<Type>();

        tmp.x = this->x / vec.x;
        tmp.y = this->y / vec.y;
        tmp.z = this->z / vec.z;

        return tmp;
    }
};

/**
 * @typedef define a typedef for 3D float vector
 */
typedef Vector3<float> Vector3f;

#endif //PROJECT_VECTOR_HPP
