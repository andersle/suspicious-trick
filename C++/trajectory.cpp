#include "trajectory.h"
#include <iostream>

/** Open the specified trajectory file.
 * \param filename The name of the file containing the trajectory.
 * \param properties List of the properties to expect for each atom. 
 * Must be in the correct order!
 */
Trajectory::Trajectory(const std::string& filename, 
		const std::vector<Atoms::Property>& properties)
	: filename(filename),
	properties(properties)
{
	file.open(filename.c_str(), std::ios::binary);

	//for the properties stored as 3-vectors, need to know which property
	//triggers us having reached the point at which we can add the 3-vector
	//this is just whichever comes last in the properties vector
	for (auto it = properties.rbegin(); it != properties.rend(); ++it) {
		switch (*it) {
			case Atoms::Property::X:
			case Atoms::Property::Y:
			case Atoms::Property::Z:
				if (ppt.x == Atoms::Property::NULL_PROPERTY) {
					ppt.x = *it;
				}
				break;
			case Atoms::Property::XS:
			case Atoms::Property::YS:
			case Atoms::Property::ZS:
				if (ppt.xs == Atoms::Property::NULL_PROPERTY) {
					ppt.xs = *it;
				}
				break;
			case Atoms::Property::XSU:
			case Atoms::Property::YSU:
			case Atoms::Property::ZSU:
				if (ppt.xsu == Atoms::Property::NULL_PROPERTY) {
					ppt.xsu = *it;
				}
				break;
			case Atoms::Property::XU:
			case Atoms::Property::YU:
			case Atoms::Property::ZU:
				if (ppt.xu == Atoms::Property::NULL_PROPERTY) {
					ppt.xu = *it;
				}
				break;
			case Atoms::Property::VX:
			case Atoms::Property::VY:
			case Atoms::Property::VZ:
				if (ppt.v == Atoms::Property::NULL_PROPERTY) {
					ppt.v = *it;
				}
				break;
			case Atoms::Property::FX:
			case Atoms::Property::FY:
			case Atoms::Property::FZ:
				if (ppt.f == Atoms::Property::NULL_PROPERTY) {
					ppt.f = *it;
				}
				break;
			case Atoms::Property::IX:
			case Atoms::Property::IY:
			case Atoms::Property::IZ:
				if (ppt.i == Atoms::Property::NULL_PROPERTY) {
					ppt.i = *it;
				}
				break;
			default:
				break;
		}	
	}
}

/** Read a single frame from the trajectory.
 * \return Atoms object populated with all  data about the timestep, including
 * an Atoms::error flag which the user must check to ensure that no errors
 * ocurred during the read.
 */
Atoms Trajectory::readFrame()
{
	Atoms a;

	if (!file.is_open()) {
		a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}

	file.read(ubi.buf, sizeof(int64_t));
	if (file.fail()) {
		if (file.eof())
			a.errorflag = Atoms::error::END_OF_FILE;
		else
			a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}
	a.timestep = ubi.i;
	
	file.read(ubi.buf, sizeof(int64_t));
	if (file.fail()) {
		a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}
	a.n = ubi.i;
	
	// Reserve enough memory for the properties that we want and the number
	// of atoms we're about to read. By reserving, we save a little time on
	// allocation later (but not much, STL allocation is quick), and also
	// move any out of memory errors to the start of the read process
	// (hopefully).
	// If the file contains too many atoms, we'll get a bad_alloc exception,
	// which we leave the caller of this function to deal with.
	for (auto p : properties) {
		switch (p) {
			case Atoms::Property::ID:
				a.id.reserve(a.n);
				break;
			case Atoms::Property::TYPE:
				a.type.reserve(a.n);
				break;
			case Atoms::Property::MOL:
				a.mol.reserve(a.n);
				break;
			case Atoms::Property::MASS:
				a.mass.reserve(a.n);
				break;
			case Atoms::Property::X:
				a.x.reserve(a.n);
				break;
			case Atoms::Property::Y:
				a.x.reserve(a.n);
				break;
			case Atoms::Property::Z:
				a.x.reserve(a.n);
				break;
			case Atoms::Property::XS:
				a.xs.reserve(a.n);
				break;
			case Atoms::Property::YS:
				a.xs.reserve(a.n);
				break;
			case Atoms::Property::ZS:
				a.xs.reserve(a.n);
				break;
			case Atoms::Property::XU:
				a.xu.reserve(a.n);
				break;
			case Atoms::Property::YU:
				a.xu.reserve(a.n);
				break;
			case Atoms::Property::ZU:
				a.xu.reserve(a.n);
				break;
			case Atoms::Property::XSU:
				a.xsu.reserve(a.n);
				break;
			case Atoms::Property::YSU:
				a.xsu.reserve(a.n);
				break;
			case Atoms::Property::ZSU:
				a.xsu.reserve(a.n);
				break;
			case Atoms::Property::IX:
				a.image_flags.reserve(a.n);
				break;
			case Atoms::Property::IY:
				a.image_flags.reserve(a.n);
				break;
			case Atoms::Property::IZ:
				a.image_flags.reserve(a.n);
				break;
			case Atoms::Property::VX:
				a.v.reserve(a.n);
				break;
			case Atoms::Property::VY:
				a.v.reserve(a.n);
				break;
			case Atoms::Property::VZ:
				a.v.reserve(a.n);
				break;
			case Atoms::Property::FX:
				a.f.reserve(a.n);
				break;
			case Atoms::Property::FY:
				a.f.reserve(a.n);
				break;
			case Atoms::Property::FZ:
				a.f.reserve(a.n);
				break;
			case Atoms::Property::Q:
				a.q.reserve(a.n);
				break;
			case Atoms::Property::NULL_PROPERTY:
				//NULL_PROPERTY needs no handling
				break;
		}
	}

	file.read(ui.buf, sizeof(int));
	if (file.fail()) {
		a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}	
	if (ui.i != 0) {
		a.errorflag = Atoms::error::TRICLINIC_BOX;
		return a;
	}

	for (int j = 0; j < 2; ++j) {
		for (int i = 0; i < 3; ++i) {
			file.read(ui.buf, sizeof(int));
			if(file.fail()) {
				a.errorflag = Atoms::error::FILE_ERROR;
				return a;
			}
			if (ui.i == 0)
				a.boxboundaries[i][j] = 'p';
			else if (ui.i == 1)
				a.boxboundaries[i][j] = 'f';
			else if (ui.i == 2)
				a.boxboundaries[i][j] = 's';
			else if (ui.i == 3)
				a.boxboundaries[i][j] = 'm';
			else {
				a.errorflag = Atoms::error::BAD_BOUNDARY;
				return a;
			}
		}
	}

	std::array<double, 6> box;
	for (int i = 0; i < 6; ++i) {
		file.read(ud.buf, sizeof(double));
		if (file.fail()) {
			a.errorflag = Atoms::error::FILE_ERROR;
			return a;
		}
		box[i] = ud.d;
	}
	a.box_lo[0] = box[0];
	a.box_lo[1] = box[2];
	a.box_lo[2] = box[4];
	a.box_hi[0] = box[1];
	a.box_hi[1] = box[3];
	a.box_hi[2] = box[5];

	file.read(ui.buf, sizeof(int));
	if (file.fail()) {
		a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}

	a.num_fields = static_cast<unsigned int>(ui.i);
	if (a.num_fields != properties.size()) {
		a.errorflag = Atoms::error::BAD_PROPERTY_COUNT;
		return a;
	}

	file.read(ui.buf, sizeof(int));
	if (file.fail()) {
		a.errorflag = Atoms::error::FILE_ERROR;
		return a;
	}
	int nprocs = ui.i; //number of processors used

	for (int i = 0; i < nprocs; ++i) {
		file.read(ui.buf, sizeof(int));
		if (file.fail()) {
			a.errorflag = Atoms::error::FILE_ERROR;
			return a;
		}
		int bufsize = ui.i; //number of doubles that follow
		if (bufsize % a.num_fields != 0) {
			// the number of atoms in this block is
			// bufsize/num_fields. if this isn't an integer,
			// something has gone badly wrong
			a.errorflag = Atoms::error::FILE_CORRUPT;
			return a;
		}
		int atoms_in_block = bufsize / a.num_fields;
		//the vector controls a contiguous memory chunk of size
		//bufsize*sizeof(double)
		std::vector<double_> buffer;
		buffer.resize(bufsize);
		file.read(buffer.data()->buf, bufsize*sizeof(double));
		//now unpack this buffer
		for (int j = 0; j < atoms_in_block; ++j) {
			Atoms::Vect3<double> x;
			Atoms::Vect3<double> xs;
			Atoms::Vect3<double> xsu;
			Atoms::Vect3<double> xu;
			Atoms::Vect3<double> v;
			Atoms::Vect3<double> f;
			Atoms::Vect3<int> i;

			for (unsigned int k = 0; k < a.num_fields; ++k) {
				double val = buffer[j*a.num_fields + k].d;
				switch (properties[k]) {
					case Atoms::Property::ID:
						a.id.emplace_back(static_cast<int>(val));
						break;
					case Atoms::Property::TYPE:
						a.type.emplace_back(static_cast<int>(val));
						break;
					case Atoms::Property::MOL:
						a.mol.emplace_back(static_cast<int>(val));
						break;
					case Atoms::Property::MASS:
						a.mass.emplace_back(val);
						break;
					case Atoms::Property::X:
						x.x = val;
						if (ppt.x == Atoms::Property::X)
							a.x.emplace_back(x);
						break;
					case Atoms::Property::Y:
						x.y = val;
						if (ppt.x == Atoms::Property::Y)
							a.x.emplace_back(x);
						break;
					case Atoms::Property::Z:
						x.z = val;
						if (ppt.x == Atoms::Property::Z)
							a.x.emplace_back(x);
						break;
					case Atoms::Property::XS:
						xs.x = val;
						if (ppt.xs == Atoms::Property::XS)
							a.xs.emplace_back(xs);
						break;
					case Atoms::Property::YS:
						xs.y = val;
						if (ppt.xs == Atoms::Property::YS)
							a.xs.emplace_back(xs);
						break;
					case Atoms::Property::ZS:
						xs.z = val;
						if (ppt.xs == Atoms::Property::ZS)
							a.xs.emplace_back(xs);
						break;
					case Atoms::Property::XSU:
						xsu.x = val;
						if (ppt.xsu == Atoms::Property::XSU)
							a.xsu.emplace_back(xsu);
						break;
					case Atoms::Property::YSU:
						xsu.y = val;
						if (ppt.xsu == Atoms::Property::YSU)
							a.xsu.emplace_back(xsu);
						break;
					case Atoms::Property::ZSU:
						xsu.z = val;
						if (ppt.xsu == Atoms::Property::ZSU)
							a.xsu.emplace_back(xsu);
						break;
					case Atoms::Property::XU:
						xu.x = val;
						if (ppt.xu == Atoms::Property::XU)
							a.xu.emplace_back(xu);
						break;
					case Atoms::Property::YU:
						xu.y = val;
						if (ppt.xu == Atoms::Property::YU)
							a.xu.emplace_back(xu);
						break;
					case Atoms::Property::ZU:
						xu.z = val;
						if (ppt.xu == Atoms::Property::ZU)
							a.xu.emplace_back(xu);
						break;
					case Atoms::Property::VX:
						v.x = val;
						if (ppt.v == Atoms::Property::VX)
							a.v.emplace_back(v);
						break;
					case Atoms::Property::VY:
						v.y = val;
						if (ppt.v == Atoms::Property::VY)
							a.v.emplace_back(v);
						break;
					case Atoms::Property::VZ:
						v.z = val;
						if (ppt.v == Atoms::Property::VZ)
							a.v.emplace_back(v);
						break;
					case Atoms::Property::FX:
						f.x = val;
						if (ppt.f == Atoms::Property::FX)
							a.f.emplace_back(f);
						break;
					case Atoms::Property::FY:
						f.y = val;
						if (ppt.f == Atoms::Property::FY)
							a.f.emplace_back(f);
						break;
					case Atoms::Property::FZ:
						f.z = val;
						if (ppt.f == Atoms::Property::FZ)
							a.f.emplace_back(f);
						break;
					case Atoms::Property::IX:
						i.x = static_cast<int>(val);
						if (ppt.i == Atoms::Property::IX)
							a.image_flags.emplace_back(i);
						break;
					case Atoms::Property::IY:
						i.y = val;
						if (ppt.i == Atoms::Property::IY)
							a.image_flags.emplace_back(i);
						break;
					case Atoms::Property::IZ:
						i.z = val;
						if (ppt.i == Atoms::Property::IZ)
							a.image_flags.emplace_back(i);
						break;
					case Atoms::Property::Q:
						a.q.emplace_back(val);
						break;
					case Atoms::Property::NULL_PROPERTY:
						//NULL_PROPERTY needs no handling
						break;
				}	
			}
		}
		
	}
	return a;
}
