#ifndef __DEBUG_LOG_H__
#define __DEBUG_LOG_H__
#include <string>
#include <iostream>
#include <memory>
#include <fstream>
#include <vector>
#include <omp.h>

enum class Log_Destination {
	CONSOLE,
	FILE
};

class Log {
public:
	static std::shared_ptr<Log> get_instance(std::string log_filename);
  static std::shared_ptr<Log> get_instance(std::string log_filename, std::string code_file);

	void set_destination(std::string log_filename,
						   Log_Destination destination);

	void Log_this(std::string code_file, int code_line, std::string message);
	void Log_this(std::string message, int code_line);


private:
	Log_Destination m_destination;
	std::ofstream m_log_stream;
  std::string m_log_filename;
  std::string m_code_file;
	static std::shared_ptr<Log> m_instance;

	void log_message(const std::string& message);

};

// extern std::vector<std::shared_ptr<Log>> debug_logs;
#endif // !__DEBUG_LOG__
