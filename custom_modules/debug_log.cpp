#include "./debug_log.h"
#include <memory>
#include <string>

std::shared_ptr<Log> Log::m_instance;
void Log::set_destination(std::string log_filename = "",
					   Log_Destination destination = Log_Destination::CONSOLE) {
	m_destination = destination;

	if (destination == Log_Destination::FILE && !log_filename.empty()) {
		m_log_stream.open(log_filename);
		if (!m_log_stream.good()) {
			std::cerr << "Can't Open Log File" << std::endl;
			m_destination = Log_Destination::CONSOLE;
      m_log_filename="";
		}
    m_log_filename=log_filename;
    m_log_stream.close();
    
	}
}

std::shared_ptr<Log> Log::get_Log(std::string log_filename) {
	if (m_instance == nullptr) {
    m_instance.reset(new Log());
	}
  m_instance->set_destination(log_filename, Log_Destination::FILE);
	return m_instance;
}
std::shared_ptr<Log> Log::get_Log(std::string log_filename, std::string code_file) {
	if (m_instance == nullptr) {
    m_instance.reset(new Log());
	}
  m_instance->set_destination(log_filename, Log_Destination::FILE);
  m_instance->m_code_file=code_file;
	return m_instance;
}

void Log::Log_this(std::string code_file, int code_line, std::string message) {
  #pragma omp critical
  {	
		code_file += " : " + std::to_string(code_line) + " : ";
		message = code_file + message;
	  log_message(message);
  }
}
void Log::Log_this(std::string message, int code_line) {
  #pragma omp critical
  {
	  std::string code_file= m_code_file;
		code_file += " : before " + std::to_string(code_line) + " : ";
		message = code_file + message;
	  log_message(message);
    // debug_logs.push_back(m_instance);
  }
}

void Log::log_message(const std::string& message) {
	if (m_destination == Log_Destination::FILE) {
    #pragma omp critical
    {
      m_log_stream.open(m_log_filename, std::ios_base::app);
      m_log_stream << message << std::endl;
      m_log_stream.close();
    }
	}
	else {
		std::cout << message << std::endl;
	}
}
