#include "header.hpp"

/* Displacement ***************************************************************/

void
Displacement::flush ()
{
	for (int j=0; j<i; j++) {
		of << m_buffer[j] << std::endl;
	}
}

void
Displacement::insert (double y)
{
	m_buffer[i] = y;
	i++;
	if (i==BUFF_SIZE) {
		this->flush();
		i=0;
	}
}

void
Displacement::close ()
{
	this->of.close();
}
