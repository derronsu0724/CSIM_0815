#include <string>  
#include <boost/regex.hpp>
#include <spdlog/spdlog.h> // NOLINT
#include <spdlog/sinks/rotating_file_sink.h>  // NOLINT
#include "spdlog/sinks/stdout_color_sinks.h"  // NOLINT
#include <filesystem>
#include "log.h"
void testSpdLog(std::string str1) {
    const std::string logFileName = "log.log"; // log文件名      
    // 检查文件是否存在，如果存在则删除  
    if (std::filesystem::exists(logFileName)) {  
        std::filesystem::remove(logFileName);  
    }  
	// 文件日志定义，设定日志最大100M，且最多保留10个
	auto fileLogger = spdlog::rotating_logger_mt("fileLogger", logFileName, 1024 * 1024 * 100, 10);
	//大于等于该等级的将被输出
	fileLogger->set_level(spdlog::level::trace);
 
	//定义输出格式
	spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]");  //log类型 文件->函数名->行数   输出内容
	spdlog::set_default_logger(fileLogger);
 
	//SPDLOG_LOGGER_TRACE(fileLogger, "trace输出：{}-{}", i, "测试trace");
	//SPDLOG_LOGGER_DEBUG(fileLogger, "debug输出：{}-{}", i, "测试debug");
	SPDLOG_LOGGER_INFO(fileLogger, "info:{}", str1);
	//SPDLOG_LOGGER_WARN(fileLogger, "warn输出：{}-{}", i, "测试warn");
	//SPDLOG_LOGGER_ERROR(fileLogger, "Some error message");
	//SPDLOG_LOGGER_CRITICAL(fileLogger, "Some critical message");

}

void testMultiLog(std::string str1) {
    const std::string logFileName = "logs.log"; // log文件名      
    // 检查文件是否存在，如果存在则删除  
    if (std::filesystem::exists(logFileName)) {  
        std::filesystem::remove(logFileName);  
    }  
	//文件sink
	auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_st>(logFileName, 1024 * 1024 * 100, 10);
	file_sink->set_level(spdlog::level::debug);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]"); 
	//控制台sink
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
	console_sink->set_level(spdlog::level::debug);
	console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]"); 
	std::vector<spdlog::sink_ptr> sinks;
	sinks.push_back(console_sink);
	sinks.push_back(file_sink);
 
	auto multiLogger = std::make_shared<spdlog::logger>("multiSink", begin(sinks), end(sinks));
	multiLogger->set_level(spdlog::level::debug);
	spdlog::set_default_logger(multiLogger);
	//SPDLOG_LOGGER_TRACE(multiLogger, "trace输出：{}-{}", i, "测试trace");
	//SPDLOG_LOGGER_DEBUG(multiLogger, "debug输出：{}-{}", i, "测试debug");
	SPDLOG_LOGGER_INFO(multiLogger, "info:{}", str1);
	SPDLOG_LOGGER_WARN(multiLogger, "warn:{}", str1);
	//SPDLOG_LOGGER_ERROR(multiLogger, "Some error message");
	//SPDLOG_LOGGER_CRITICAL(multiLogger, "Some critical message");

}
