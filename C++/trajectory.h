#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <cstdint>
#include <fstream>
#include <string>

#include "atoms.h"

/** Reads data from trajectory files.
 * Trajectory reads a LAMMPS dump file one step at a time, returning an Atoms
 * object after each frame. It supports only binary dump files at the moment.
 * The user must supply a list of all the fields to expect from the dump file,
 * and the order in which to expect them.
 */
class Trajectory {
public:

	Trajectory(const std::string&, const std::vector<Atoms::Property>&);

	Atoms readFrame(); 
private:
	/** Used for reading the LAMMPS bigint type */
	union bigint_ {
		char buf[sizeof(int64_t)];
		int64_t i;
	} ubi;

	/** Used to read a standard double */
	union double_ {
		char buf[sizeof(double)];
		double d;
	} ud;

	/** Used to read a standard int */
	union int_ {
		char buf[sizeof(int)];
		int i;
	} ui;

	/** Properties which cause a Vect3 struct to be pushed onto the relevant
	 * vector. */
	struct PropertyPushTriggers { 
		Atoms::Property x;
		Atoms::Property xs;
		Atoms::Property xsu;
		Atoms::Property xu;
		Atoms::Property v;
		Atoms::Property f;
		Atoms::Property i;

		PropertyPushTriggers()
			: x(Atoms::Property::NULL_PROPERTY),
			xs(Atoms::Property::NULL_PROPERTY),
			xsu(Atoms::Property::NULL_PROPERTY),
			xu(Atoms::Property::NULL_PROPERTY),
			v(Atoms::Property::NULL_PROPERTY),
			f(Atoms::Property::NULL_PROPERTY),
			i(Atoms::Property::NULL_PROPERTY) {}
	} ppt;

	/** Name of the trajectory file to read from. */
	const std::string filename;
	/** List of the properties to read for each atom */
	const std::vector<Atoms::Property>& properties;
	std::ifstream file;
};

#endif
