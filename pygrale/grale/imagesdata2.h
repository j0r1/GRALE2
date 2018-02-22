#ifndef IMAGESDATA2_H

#define IMAGESDATA2_H

#include <grale/imagesdata.h>

class ImagesData2 : public grale::ImagesData
{
public:
	bool read2(serut::SerializationInterface &si) { return ImagesData::read(si); }
	bool write2(serut::SerializationInterface &si) const { return ImagesData::write(si); }
};

#endif // IMAGESDATA2_H
