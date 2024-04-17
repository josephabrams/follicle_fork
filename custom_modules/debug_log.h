#ifndef __DEBUG_LOG_H__
#define __DEBUG_LOG_H__
#include <string>
#include <iostream>
#include <memory>
#include <fstream>


enum class Log_Destination {
	CONSOLE,
	FILE
};

class Log {
public:
	static std::shared_ptr<Log> get_instance();

	void set_destination(std::string log_filename,
						   Log_Destination destination);

	void Log_this(std::string code_file, int code_line, std::string message);
	void Log_this(std::string message);


private:
	Log_Destination m_destination;
	std::ofstream m_log_stream;
  std::string m_log_filename;
	static std::shared_ptr<Log> m_instance;

	void log_message(const std::string& message);

};

#endif // !__DEBUG_LOG__
