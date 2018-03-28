#ifndef COMMONGA_H

#define COMMONGA_H

template<class GABase>
class CommonGA : public GABase
{
	void feedbackStatus(const std::string &str) const override
	{
		//cerr << "cerr:GAFEEDBACK:" << str << endl;
		WriteLineStdout("GAFEEDBACK:" + str);
	}

	void onMessage(const std::string &msg) override
	{
		//cerr << "cerr:GAMESSAGESTR:" << msg << endl;
		WriteLineStdout("GAMESSAGESTR:" + msg);
	}

	void onMessage(const std::vector<uint8_t> &msg) override
	{
		stringstream ss;

		ss << "GAMESSAGEBYTES:" << msg.size();

		//cerr << "cerr:" << ss.str() << endl;
		WriteLineStdout(ss.str());
		if (msg.size() > 0)
			WriteBytesStdout(msg);
	}
};

#endif // COMMONGA_H
