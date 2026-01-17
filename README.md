# 3Dmol - 3D分子结构可视化工具

使用 RDKit 和 py3Dmol 库生成和可视化 3D 分子结构的工具。

## 项目概述

本项目提供了一个 Python 脚本，可以自动生成分子的 3D 结构模型，并将其保存为交互式 HTML 文件。当前示例生成 H₂S（硫化氢）分子的 3D 结构。

## 文件说明

### 3Dmolecule.py
主程序文件，用于生成中性分子的 3D 结构，包含以下功能：
- 使用 SMILES 标记法定义分子结构
- 自动添加隐含氢原子
- 通过 ETKDG 算法生成 3D 坐标
- 使用 MMFF 力场优化分子几何结构
- 生成交互式 3D 可视化 HTML 文件

**输出文件**：根据 SMILES 生成相应的 HTML 文件（例如 `(CH3)2S_molecule.html`）

### 3D-ion.py
离子结构生成脚本，用于生成离子的 3D 结构，功能包括：
- 使用带电荷的 SMILES 标记法定义离子结构
- 支持带正电和负电的离子
- 生成离子的 3D 几何结构
- 生成交互式 3D 可视化 HTML 文件

**示例**：H-O-Zn-O 离子（`O[Zn+2][O-]`）生成 `H_O_Zn_O_ion.html`

### 3D-adsorption.py
吸附结构生成脚本，用于生成分子或离子之间的吸附结构，功能包括：
- 同时加载两个或多个分子/离子
- 设置合理的空间位置以显示吸附配置
- 支持不同类型的吸附相互作用：
  - S 原子与 Zn 原子的配位吸附
  - O 原子与 C 原子的相互作用
- 生成显示吸附结构的交互式 3D 可视化 HTML 文件

**示例**：H-O-Zn-O 离子与 CH3SH 分子的吸附结构，生成 `H_O_Zn_O_ion_CH3SH_adsorption.html`

### HTML 输出文件
所有脚本生成的 HTML 文件都包含可交互的 3D 分子/离子/吸附结构可视化，可直接在网页浏览器中打开。

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

### 1. 生成中性分子结构

运行脚本生成中性分子的 3D 结构：

```bash
python3 3Dmolecule.py
```

脚本会根据 SMILES 字符串生成对应的 HTML 文件（例如 `(CH3)2S_molecule.html`）。

### 2. 生成离子结构

运行脚本生成离子的 3D 结构：

```bash
python3 3D-ion.py
```

脚本会生成包含离子结构的 HTML 文件（例如 `H_O_Zn_O_ion.html`）。

**修改离子结构**：编辑 `3D-ion.py` 中的 SMILES 字符串，使用带电荷的格式：

```python
smiles_string = 'O[Zn+2][O-]'  # Zn²⁺ 离子，一个 O⁻ 负离子
```

### 3. 生成吸附结构

运行脚本生成分子与离子的吸附结构模型：

```bash
python3 3D-adsorption.py
```

脚本会生成展示吸附配置的 HTML 文件（例如 `H_O_Zn_O_ion_CH3SH_adsorption.html`）。

**修改吸附结构**：编辑 `3D-adsorption.py` 中的 SMILES 字符串来改变吸附分子和吸附剂离子。

### 查看结果

使用任何现代网页浏览器打开生成的 HTML 文件：

```bash
# 在浏览器中打开
open H_O_Zn_O_ion.html              # macOS
xdg-open H_O_Zn_O_ion.html          # Linux
start H_O_Zn_O_ion.html             # Windows
```

或在文件管理器中双击打开。

## 修改分子/离子/吸附结构

### 修改中性分子结构

编辑 `3Dmolecule.py` 中的 SMILES 字符串：

```python
smiles_string = '(CH3)2S'  # 二甲硫醚
```

### 修改离子结构

编辑 `3D-ion.py` 中的 SMILES 字符串，包括离子的电荷信息：

```python
smiles_string = 'O[Zn+2][O-]'  # H-O-Zn-O 离子
```

### 修改吸附结构

编辑 `3D-adsorption.py` 中的分子和离子定义：

```python
smiles_CH3SH = 'CS'           # 被吸附的分子（甲硫醇）
smiles_ion = 'O[Zn+2][O-]'    # 吸附剂（H-O-Zn-O 离子）
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

### 中性分子结构生成（3Dmolecule.py）
1. **分子定义**：使用 SMILES 表示法定义分子
2. **分子对象创建**：RDKit 解析 SMILES 并创建分子对象
3. **添加氢原子**：自动添加隐含的氢原子
4. **3D 坐标生成**：使用 ETKDG 算法生成初始 3D 坐标
5. **几何优化**：使用 MMFF 力场优化分子几何结构
6. **可视化**：使用 py3Dmol 创建交互式 3D 可视化
7. **导出**：将结果保存为 HTML 文件

### 离子结构生成（3D-ion.py）
1. **离子定义**：使用带电荷信息的 SMILES 表示法定义离子
2. **离子对象创建**：RDKit 解析 SMILES 并创建离子对象
3. **添加氢原子**：根据离子的价态自动添加氢原子
4. **3D 坐标生成**：使用 ETKDG 算法生成初始 3D 坐标
5. **几何优化**：使用 MMFF 力场优化几何结构（部分离子可能不支持）
6. **可视化**：使用自定义 HTML 生成器创建交互式 3D 可视化
7. **导出**：将结果保存为 HTML 文件

### 吸附结构生成（3D-adsorption.py）
1. **分子和离子定义**：分别定义被吸附分子和吸附剂离子
2. **独立 3D 结构生成**：为两个物种分别生成 3D 结构
3. **空间定位**：调整被吸附分子的位置以展示吸附配置
4. **结构组合**：将两个物种的坐标组合在一起
5. **可视化**：创建展示吸附相互作用的 3D 可视化
6. **导出**：将吸附结构保存为 HTML 文件

## 可视化特性

生成的 HTML 文件包含以下交互特性：
- **旋转**：鼠标左键拖动旋转分子
- **缩放**：鼠标滚轮或触摸板缩放
- **平移**：鼠标右键拖动平移视图
- **样式**：显示原子为球体，键为棍状结构

## 离子结构特性

使用 `3D-ion.py` 生成的离子结构包括：
- **电荷表示**：离子的总电荷信息会在 MOL 格式中标注
- **金属配位体**：支持过渡金属离子及其配位体
- **多原子离子**：支持复杂的多原子离子结构
- **几何表示**：清晰展示离子周围的配位环境

## 吸附结构特性

使用 `3D-adsorption.py` 生成的吸附结构包括：
- **分子对接**：展示两个物种在空间中的相对位置
- **吸附位点**：清晰标识潜在的吸附相互作用位点
- **原子类型识别**：自动识别和标记参与吸附的原子（如 S-Zn 配位、O-C 相互作用）
- **分层显示**：被吸附分子和吸附剂可单独控制的可视化表现

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
A: 编辑相应脚本中的 `setStyle()` 方法。例如：
```python
view.setStyle({}, {'sphere': {'radius': 0.4}, 'stick': {'radius': 0.2}})
```

**Q: SMILES 字符串从哪里获得？**
A: 可以从以下来源获得：
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [ChemSpider](https://www.chemspider.com/)
- [SMILES 在线生成工具](https://smiles.detritus.org/)

**Q: 如何定义带电荷的离子？**
A: 在 SMILES 字符串中使用 `[Atom+n]` 或 `[Atom-n]` 的格式：
- `[Zn+2]` 表示 Zn²⁺
- `[O-]` 表示 O⁻
- `[Na+]` 表示 Na⁺

**Q: 可以生成蛋白质或大分子结构吗？**
A: 本脚本主要用于小分子和离子。对于蛋白质等大分子，建议使用专门的工具如 PyMOL。

**Q: 吸附结构中能否显示吸附键？**
A: 当前脚本显示的是吸附分子和吸附剂在空间中的相对位置，展示潜在的相互作用位点。如需显示完整的化学键，可修改 `3D-adsorption.py` 中的分子合并代码。

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request！

## 相关资源

- [RDKit 官方文档](https://www.rdkit.org/docs/)
- [py3Dmol 官方文档](https://3dmol.csb.pitt.edu/index.html)
- [SMILES 标记法](https://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_system)
