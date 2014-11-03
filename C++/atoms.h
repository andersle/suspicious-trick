#ifndef ATOMS_H
#define ATOMS_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

/** Contains all data read from the trajectory in a given timestep. */
class Atoms {
public:
	Atoms();
	/** Various error types that might occur. */
	enum class error {
		NO_ERROR, /**< No error ocurred */
		END_OF_FILE, /**< End of file reached. */
		FILE_ERROR, /**< File not opened, unexpected EOF, etc. */
		TRICLINIC_BOX, /**< Triclinic boxes are not supported */
		BAD_BOUNDARY, /**< Unrecognised boundary type (not p,f,s,m) */
		BAD_PROPERTY_COUNT, /**< The number of properties specified by
				      the user is different to the number in the
				      datafile. */
		FILE_CORRUPT, /**< The reported buffer size for a given
				processor block isn't compatible with the 
				reported number of fields. */
	};
	/** Contains Atoms::error::NO_ERROR if nothing went wrong */
	error errorflag;
	
	/** Supported properties */
	enum class Property { 
		NULL_PROPERTY, /**< Internal placeholder value */
		ID, /**< Atom ID */
		TYPE, /**< Atom type */
		MOL, /**< Molecule to which the atom belongs */
		MASS, /**< Atom mass */
		X, /**< x coordinate */
		Y, /**< y coordinate */
		Z, /**< z coordinate */
		XS, /**< x coordinate (scaled) */
		YS, /**< y coordinate (scaled) */
		ZS, /**< z coordinate (scaled) */
		XU, /**< x coordinate (unwrapped) */
		YU, /**< y coordinate (unwrapped) */
		ZU, /**< z coordinate (unwrapped) */
		XSU, /**< x coordinate (scaled, unwrapped) */
		YSU, /**< y coordinate (scaled, unwrapped) */
		ZSU, /**< z coordinate (scaled, unwrapped) */
		IX, /**< Image flag (x direction) */
		IY, /**< Image flag (y direction) */
		IZ, /**< Image flag (z direction) */
		VX, /**< Velocity x component */
		VY, /**< Velocity y component */
		VZ, /**< Velocity z component */
		FX, /**< Force x component */
		FY, /**< Force y component */
		FZ, /**< Force z component */
		Q /**< Charge */
	};

	/** Stores 3-vectors like position, velocity and force */
	template<typename T>
	struct Vect3 {
		T x; /**< x component */
		T y; /**< y component */
		T z; /**< z component */

		inline Vect3<T> operator-(const Vect3<T>& rhs) const {
			Vect3<T> ret;
			ret.x = x - rhs.x;
			ret.y = y - rhs.y;
			ret.z = z - rhs.z;
			return ret;
		}

		inline Vect3<T> operator*(const T& factor) const {
			Vect3<T> ret;
			ret.x = x * factor;
			ret.y = y * factor;
			ret.z = z * factor;
			return ret;
		}


		inline Vect3<T> operator*(const Vect3<T>& rhs) const {
			Vect3<T> ret;
			ret.x = x * rhs.x;
			ret.y = y * rhs.y;
			ret.z = z * rhs.z;
			return ret;
		}

		inline Vect3<T> operator*=(const T& factor) {
			x *= factor;
			y *= factor;
			z *= factor;
			return *this;
		}

		inline Vect3<T> operator/(const Vect3<T>& rhs) const {
			Vect3<T> ret;
			ret.x = x / rhs.x;
			ret.y = y / rhs.y;
			ret.z = z / rhs.z;
			return ret;
		}

		inline Vect3<T> operator/(const T& rhs) const {
			Vect3<T> ret;
			ret.x = x/rhs;
			ret.y = y/rhs;
			ret.z = z/rhs;
			return ret;
		}

		inline Vect3<T> operator+(const Vect3<T>& rhs) const {
			Vect3<T> ret;
			ret.x = x + rhs.x;
			ret.y = y + rhs.y;
			ret.z = z + rhs.z;
			return ret;
		}

		inline Vect3<T> operator+=(const Vect3<T>& rhs) {
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			return *this;
		}

		inline Vect3<T> operator-=(const Vect3<T>& rhs) {
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			return *this;
		}
		inline const T len2() const {
			return x*x+y*y+z*z;
		}

		inline Vect3<T> norm() {
			T len = sqrt(len2());
			x /= len;
			y /= len;
			z /= len;

			return *this;
		}

		inline const T dot(const Vect3<T>& rhs) const {
			return x*rhs.x+y*rhs.y+z*rhs.z;
		}

		inline Vect3<T> cross(const Vect3<T>& rhs) const {
			Vect3<T> ret;
			ret.x = y*rhs.z - z*rhs.y;
			ret.y = x*rhs.z - z*rhs.x;
			ret.z = x*rhs.y - y*rhs.x;
			return ret;
		}
	};
	
	// Frame header fields

	/** Number of atoms in this container */
	uint64_t n;
	/** The timestep from which this data was taken */
	uint64_t timestep;
	/** End points of the three box axes */
	std::array<double, 3> box_hi;
	/** Start points of the three box axes */
	std::array<double, 3> box_lo;
	/** Types of the box faces: (p)eriodic, (f)ixed, ... u if unset. */
	std::array<std::array<char, 2>, 3> boxboundaries;
	/** The number of fields per atom recorded. */
	unsigned int num_fields;

	// Atom data lists
	
	/** List of atomic forces */
	std::vector<Vect3<double>> f;
	/** List of atom IDs */
	std::vector<int> id;
	/** List of atomic image flags */
	std::vector<Vect3<int>> image_flags;
	/** List of masses */
	std::vector<double> mass;
	/** List of molecule IDs */
	std::vector<int> mol;
	/** List of atomic charges */
	std::vector<double> q;
	/** List of atom types */
	std::vector<int> type;
	/** List of atomic velocities */
	std::vector<Vect3<double>> v;
	/** List of atomic positions */
	std::vector<Vect3<double>> x;
	/** List of atomic positions (scaled) */
	std::vector<Vect3<double>> xs;
	/** List of atomic positions (scaled and unwrapped) */
	std::vector<Vect3<double>> xsu;
	/** List of atomic positions (unwrapped) */
	std::vector<Vect3<double>> xu;
};

#endif
