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

std::shared_ptr<Log> Log::get_instance() {
	if (m_instance == nullptr) {
    m_instance.reset(new Log());
	}
	return m_instance;
}

void Log::Log_this(std::string code_file, int code_line, std::string message) {
	
		code_file += " : " + std::to_string(code_line) + " : ";
		message = code_file + message;
	  log_message(message);
}
void Log::Log_this(std::string message) {
	  std::string code_file= static_cast<std::string>(__FILE__);
	  std::string code_line= std::to_string(__LINE__);
		code_file += " : before " + code_line + " : ";
		message = code_file + message;
	  log_message(message);
}

void Log::log_message(const std::string& message) {
	if (m_destination == Log_Destination::FILE) {
    m_log_stream.open(m_log_filename, std::ios_base::app);
		m_log_stream << message << std::endl;
    m_log_stream.close();
	}
	else {
		std::cout << message << std::endl;
	}
}
