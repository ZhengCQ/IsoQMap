import os
import sys
import logging
from typing import Union, Optional, List

def check_file_exists(
    file_path: Union[str, List[str]],
    file_description: str = "file",
    check_readable: bool = True,
    logger: Optional[logging.Logger] = None,
    exit_on_error: bool = True
) -> bool:
    """
    通用文件存在性检查函数
    
    Args:
        file_path: 文件路径(可以是字符串或路径列表)
        file_description: 文件描述(用于错误信息)
        check_readable: 是否检查可读性
        logger: 可选的logger对象
        exit_on_error: 是否在错误时退出程序
    
    Returns:
        bool: 文件是否存在(且可读)
    
    Raises:
        FileNotFoundError: 当文件不存在且exit_on_error为False时
    """
    def _log_error(msg):
        if logger:
            logger.error(msg)
        print(f"\033[91mERROR:\033[0m {msg}", file=sys.stderr)
    
    # 处理路径列表情况
    if isinstance(file_path, list):
        missing_files = [fp for fp in file_path if not os.path.exists(fp)]
        if missing_files:
            _log_error(f"Missing {file_description} files: {', '.join(missing_files)}")
            if exit_on_error:
                sys.exit(1)
            return False
        return True
    
    # 检查单个文件
    if not os.path.exists(file_path):
        _log_error(f"{file_description.capitalize()} not found: {file_path}")
        if exit_on_error:
            sys.exit(1)
        return False
    
    if check_readable and not os.access(file_path, os.R_OK):
        _log_error(f"{file_description.capitalize()} is not readable: {file_path}")
        if exit_on_error:
            sys.exit(1)
        return False
    
    if logger:
        logger.info(f"Found {file_description}: {file_path}")
    
    return True

import logging

def setup_logger(log_file, verbose):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # 防止重复添加 handler
    if not logger.handlers:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        ch.setLevel(logging.DEBUG if verbose else logging.INFO)
        logger.addHandler(ch)