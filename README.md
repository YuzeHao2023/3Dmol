# 3Dmol - 3D分子结构可视化工具

使用 RDKit 和 py3Dmol 库生成和可视化 3D 分子结构的工具。

## 项目概述

本项目提供了一个 Python 脚本，可以自动生成分子的 3D 结构模型，并将其保存为交互式 HTML 文件。当前示例生成 H₂S（硫化氢）分子的 3D 结构。

## 文件说明

### 3Dmolecule.py
主程序文件，包含以下功能：
- 使用 SMILES 标记法定义分子结构
- 自动添加隐含氢原子
- 通过 ETKDG 算法生成 3D 坐标
- 使用 MMFF 力场优化分子几何结构
- 生成交互式 3D 可视化 HTML 文件

### H2S_molecule.html
生成的 H₂S 分子的 3D 结构可视化文件，可以在网页浏览器中打开查看。

## 系统要求

- Python 3.6+
- RDKit
- py3Dmol

## 安装依赖

```bash
pip install rdkit-pypi py3Dmol
```

或使用 conda：

```bash
conda install -c conda-forge rdkit py3Dmol
```

## 使用方法

### 运行脚本生成 H₂S 分子结构

```bash
python3 3Dmolecule.py
```

脚本会生成 `H2S_molecule.html` 文件。

### 查看结果

使用任何现代网页浏览器打开生成的 HTML 文件：

```bash
# 在浏览器中打开
open H2S_molecule.html        # macOS
xdg-open H2S_molecule.html    # Linux
start H2S_molecule.html       # Windows
```

或在文件管理器中双击打开。

## 修改分子结构

要生成其他分子的 3D 结构，编辑 `3Dmolecule.py` 中的 SMILES 字符串：

```python
smiles_string = 'S'  # H2S 硫化氢
```

### 常见分子的 SMILES 表示

| 分子 | SMILES |
|------|--------|
| 甲烷 (CH₄) | C |
| 乙烷 (C₂H₆) | CC |
| 水 (H₂O) | O |
| 二氧化碳 (CO₂) | O=C=O |
| 甲醇 (CH₃OH) | CO |
| 苯 (C₆H₆) | c1ccccc1 |

## 程序流程

1. **分子定义**：使用 SMILES 表示法定义分子
2. **分子对象创建**：RDKit 解析 SMILES 并创建分子对象
3. **添加氢原子**：自动添加隐含的氢原子
4. **3D 坐标生成**：使用 ETKDG 算法生成初始 3D 坐标
5. **几何优化**：使用 MMFF 力场优化分子几何结构
6. **可视化**：使用 py3Dmol 创建交互式 3D 可视化
7. **导出**：将结果保存为 HTML 文件

## 可视化特性

生成的 HTML 文件包含以下交互特性：
- **旋转**：鼠标左键拖动旋转分子
- **缩放**：鼠标滚轮或触摸板缩放
- **平移**：鼠标右键拖动平移视图
- **样式**：显示原子为球体，键为棍状结构

## 代码示例

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# 创建分子
smiles_string = 'S'  # H2S
mol = Chem.MolFromSmiles(smiles_string)

# 添加氢原子
mol = Chem.AddHs(mol)

# 生成 3D 结构
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)

# 创建可视化
view = py3Dmol.view(width=400, height=400)
view.addModel(Chem.MolToMolBlock(mol), 'mol')
view.setStyle({'sphere':{'radius':0.3}, 'stick':{'radius':0.1}})
view.zoomTo()

# 导出为 HTML
html_content = view.write_html()
with open('molecule.html', 'w') as f:
    f.write(html_content)
```

## 常见问题

**Q: 如何修改可视化的样式？**
A: 编辑 `3Dmolecule.py` 中的 `setStyle()` 方法。例如：
```python
view.setStyle({'cartoon': {'color': 'spectrum'}})
```

**Q: SMILES 字符串从哪里获得？**
A: 可以从以下来源获得：
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [ChemSpider](https://www.chemspider.com/)
- [SMILES 在线生成工具](https://smiles.detritus.org/)

**Q: 可以生成蛋白质或大分子结构吗？**
A: 本脚本主要用于小分子。对于蛋白质等大分子，建议使用专门的工具如 PyMOL。

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request！

## 相关资源

- [RDKit 官方文档](https://www.rdkit.org/docs/)
- [py3Dmol 官方文档](https://3dmol.csb.pitt.edu/index.html)
- [SMILES 标记法](https://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_system)
