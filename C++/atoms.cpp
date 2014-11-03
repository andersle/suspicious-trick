#include "atoms.h"

Atoms::Atoms()
	: errorflag{error::NO_ERROR},
	n{0},
	timestep{0},
	box_hi{{0.0, 0.0, 0.0}},
	box_lo{{0.0, 0.0, 0.0}},
	boxboundaries {{{{'u', 'u'}}, {{'u', 'u'}}, {{'u', 'u'}}}},
	num_fields{0}
{}
