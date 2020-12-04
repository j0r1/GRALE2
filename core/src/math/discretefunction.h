#pragma once

#include "graleconfig.h"
#include <errut/booltype.h>
#include <vector>

namespace grale
{

template <class T>
class DiscreteFunction
{
public:
	DiscreteFunction()
	{
		m_values = { (T)0, (T)0 };
		m_x0 = (T)0;
		m_x1 = (T)1;
		m_diff = (T)1;
	}
	~DiscreteFunction() { }
	errut::bool_t init(T x0, T x1, const std::vector<T> &values);

	T operator()(T x) const;

	const std::vector<T> &getValues() const { return m_values; }
	void getLimits(T &x0, T &x1) const { x0 = m_x0; x1 = m_x1; }
private:
	std::vector<T> m_values;
	T m_x0, m_x1, m_diff;
};

template <class T>
inline errut::bool_t DiscreteFunction<T>::init(T x0, T x1, const std::vector<T> &values)
{
	if (x0 >= x1)
		return "x0 must be smaller than x1";
	if (values.size() < 2)
		return "At least two values must be present";
	
	m_values = values;
	m_x0 = x0;
	m_x1 = x1;
	m_diff = (m_x1 - m_x0)/((T)(m_values.size()-1));
	return true;
}

template <class T>
inline T DiscreteFunction<T>::operator()(T x) const
{
	T fi = (x - m_x0)/m_diff;
	if (fi < 0)
		return m_values[0];
	int i = (int)fi;
	if (i >= (int)m_values.size()-1)
		return m_values[m_values.size()-1];
	
	T di = fi - (T)i;
	return m_values[i] * ((T)1 - di) + m_values[i+1]*di;
}

}
