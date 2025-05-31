import os
import sys
from pathlib import Path
from typing import List, Optional

class BinPathFinder:
    """找到软件主目录"""
    def __init__(self, package_name: str):
        self.package_name = package_name
        self.search_paths = self._init_search_paths()
    
    def _init_search_paths(self) -> List[Path]:
        """初始化搜索路径（包含Conda环境支持）"""
        paths = []
        
        # 1. 开发环境路径（最高优先级）
        dev_path = Path(__file__).absolute().parent.parent
        if dev_path.exists():
            paths.append(dev_path)
        
        
        # 2. Conda环境路径（如果检测到Conda）
        if 'CONDA_PREFIX' in os.environ:
            conda_paths = [
                Path(os.environ['CONDA_PREFIX']) / 'share' / self.package_name,
                Path(os.environ['CONDA_PREFIX']) / 'lib' / self.package_name
            ]
            paths.extend(p for p in conda_paths if p.exists())
        
        # 3. 用户目录路径
        user_path = Path.home() / '.local' / 'share' / self.package_name
        if user_path.exists():
            paths.append(user_path)
        
        # 4. 系统级路径
        sys_paths = [
            Path('/usr/local/share') / self.package_name,
            Path('/usr/share') / self.package_name
        ]
        paths.extend(p for p in sys_paths if p.exists())
        return paths
    
    def find(self, relative_path: str) -> Optional[Path]:
        """查找相对路径文件"""
        for base_path in self.search_paths:
            target = base_path / relative_path
            if target.exists():
                return target
        return None

    def find_all(self, relative_path: str) -> List[Path]:
        """查找所有匹配路径（按优先级排序）"""
        found = []
        for base_path in self.search_paths:
            target = base_path / relative_path
            if target.exists():
                found.append(target)
        return found
