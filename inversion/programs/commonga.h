#ifndef COMMONGA_H

#define COMMONGA_H

template<class GABase>
class CommonGA : public GABase
{
	void feedbackStatus(const std::string &str) const override
	{
		cerr << "GA feedback: " << str << endl;
	}
	void onCurrentBest(const std::list<mogal::Genome *> &bestgenomes) const override
	{
		cerr << "Current best: ";
		for (auto g : bestgenomes)
			cerr << "( " << g->getFitnessDescription() << ")";
		cerr << endl;
	}
};

#endif // COMMONGA_H
