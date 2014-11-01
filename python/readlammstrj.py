"""
A simple class to read a binary lamms file
"""
import struct

class lammps_trj(object):
    def __init__(self, filename):
        self.variables_order = ['timestep', 'ndump', 'triclinic',
                                'xlo-flag', 'xhi-flag', 'ylo-flag', 'yhi-flag', 
                                'zlo-flag', 'zhi-flag', 'xlo', 'xhi', 
                                'yho', 'yhi','zlo','zhi','fields per atom',
                                'number of blocks']
        self.head_fmt = '<qqiiiiiiiddddddii' # format for reading the header
        self.head_buff = 100 # buffer length for the header
        self.filename = filename
        self.int_size = 4
        self.double_size = 8

    def read_frame(self):
        """
        Function to read a frame from a binary lammps file
        """
        data = {}
        this_chunk = 0
        buff = self.f.read(self.head_buff)
        if not buff: return data, this_chunk
        x = struct.unpack(self.head_fmt, buff)
        for bi,v in zip(x, self.variables_order):
            data[v] = bi
        this_chunk += self.head_buff
        data['blockdata'] = []
        for x in range(data['number of blocks']): # one block per cpu
            buff_size = struct.unpack('i', self.f.read(self.int_size))[0] # total number of doubles to follow
            this_chunk += self.int_size
            fmt = 'd'*buff_size
            thisblock = struct.unpack(fmt, self.f.read(buff_size*self.double_size))
            this_chunk += (buff_size*self.double_size)
            data['blockdata'].append(thisblock)
        return data, this_chunk

    def read_trajectory(self):
        self.f = open(self.filename, 'rb')
        while True:
            data, new_chunk, = self.read_frame()
            if not data:
                break
            yield data
        self.f.close()

if __name__=='__main__':
    lt = lammps_trj('track.lammpstrj.bin')
    for frame in lt.read_trajectory():
        print 'Step: {}'.format(frame['timestep'])
        #for b in frame['blockdata']: 
        #    m = np.reshape(b,(-1,5))
        #    x,y,z = m[:,2],m[:,3],m[:,4]
        

