# Codex 云端 Linux MATLAB

当前交付状态：安装与运行脚本已准备；本次会话缺失的四个运行库已下载并解包到用户目录。MATLAB 尚未安装，学校许可接入及完整分析尚未实际验证。

本环境为 Debian 13、x86_64、2 核、8 GiB 内存，采用原生 MATLAB R2026b。安装路径默认是 `/workspace/MATLAB/R2026b`，可用 `MATLAB_ROOT` 指定其他路径。缺失系统库由 `runtime-deps.sh` 在无 root 权限下准备到 `/workspace/MATLAB/system-libs`，可用 `MATLAB_SYSTEM_LIBS` 修改位置。

## 安装前需要完成

- 在云环境的网络设置中允许 MathWorks 下载域名。当前策略未放行，访问会返回 403；至少需要覆盖 `www.mathworks.com`、`ssd.mathworks.com` 及下载过程中使用的其他 MathWorks 子域名。修改本地网络策略文件不能放行网络，也不要绕过会话代理。
- 取得学校网络许可服务器信息：`lmgrd` 地址和端口，以及 MATLAB vendor daemon 的端口。两者都需要云环境 TCP 访问许可；若服务器只允许校内访问，还需要已配置的学校 VPN 或其他获准网络连接。当前会话未配置 VPN，TCP 网关也无法连接。学校地址仅存入会话工作目录，不提交到 GitHub。
- MATLAB 网络许可证不支持普通 HTTP/HTTPS 代理；设置下载代理不能使许可证连接自动可用。本平台的私网 TCP 接入还需平台 VPN、服务器域名 TCP 授权及与许可证协议兼容的转发，由运维或学校管理员确认；本平台没有可直接从执行器访问任意公网 TCP 端口的通道。
- 确认学校许可证允许当前 MATLAB 版本和云端使用。将服务器端点配置到环境变量 `MLM_LICENSE_FILE`；不要把账号密码或许可文件内容提交到仓库。
- 当前 `/workspace` 约有 7.5 GiB、独立 `/tmp` 约有 4.9 GiB 可用。仅安装 MATLAB 本体、保留默认 JRE，临时文件放 `/tmp`。安装峰值空间是否足够仍须实际验证；需要其他工具箱时应扩容。

## 安装和运行

以下命令在仓库根目录执行。先用学校提供的实际端点替换示例：

```bash
export MLM_LICENSE_FILE='27000@license.example.edu'
bash cloud-matlab/bootstrap.sh --check
bash cloud-matlab/bootstrap.sh
bash cloud-matlab/matlab.sh -batch "disp(version)"
bash cloud-matlab/matlab.sh -batch "run('cloud-matlab/run_analysis.m')"
```

`--check` 用于检查安装条件，不等于已经安装或授权成功。`runtime-deps.sh` 负责缺失运行库；`matlab.sh` 使用上述 MATLAB 和系统库路径启动程序。整个流程不要求 `sudo`，保留环境提供的代理和证书设置。

本会话的默认用户目录不可写，启动脚本将 MATLAB 偏好目录设置为 `/workspace/work/matlab-preferences/R2026b`；可通过 `MATLAB_PREFDIR` 覆盖，不修改 `HOME`。

在本次会话，提供的许可端点保存在 `/workspace/work/matlab-bootstrap/license-server.txt`（权限 600）；没有显式设置 `MLM_LICENSE_FILE` 时，两个启动脚本会读取该文件。新环境需要重新设置端点或通过云环境的运行变量传入。端点存在不代表服务器已可达或许可已经验证。

缺失依赖从 Debian 官方签名仓库下载，APT 索引写入工作目录，不修改系统 APT 状态。脚本保留系统的 CA 校验和下载代理。使用前可执行 `bash -n cloud-matlab/*.sh` 检查 shell 语法。

`run_analysis.m` 将仓库中的 `genggai1101.m` 原样复制为 `/workspace/work/matlab-validation/fdm_paper_analysis.m`，使文件名与主函数一致，仓库源码不变。运行使用当前会话搜索路径，不调用 `savepath`。它会生成图片但不显示窗口，输出位于 `/workspace/work/matlab-validation/results`。

运行成功后，脚本检查：组号 1–28 全部保留；10CF 和 15CF 各 14 组；每种材料的数据、实测排名和权重敏感性都包含相同的 14 个组号；排除组号为空；每项响应的诊断样本数为 14；汇总只有两行 `All` 情景。回归、ANOVA、嵌套 LOO、权重和 Pareto 方法不作修改。

## 官方依据

- [R2026b Linux 系统要求](https://www.mathworks.com/support/requirements/matlab-linux.html)
- [获取 MATLAB Package Manager](https://www.mathworks.com/help/install/ug/get-mpm-os-command-line.html)
- [MPM 安装参数、TMPDIR 和代理支持](https://www.mathworks.com/help/install/ug/mpminstall.html)
- [Debian 13 官方运行库清单](https://github.com/mathworks-ref-arch/container-images/blob/main/matlab-deps/r2026b/debian13/base-dependencies-amd64.txt)
- [网络许可证配置](https://www.mathworks.com/help/install/ug/network-license-files.html)
- [网络许可证不支持普通代理](https://www.mathworks.com/matlabcentral/answers/2181637-does-the-network-license-manager-support-proxy-connections)
- [网络许可证的两个 TCP 端口](https://www.mathworks.com/matlabcentral/answers/96756-what-ports-does-the-network-license-manager-use-and-how-can-i-set-those-ports)
- [MATLAB 偏好目录设置](https://www.mathworks.com/matlabcentral/answers/93696-how-do-i-change-the-matlab-preferences-directory-location)
- [MathWorks 网络域名和端口](https://www.mathworks.com/matlabcentral/answers/97265-which-ports-and-addresses-do-mathworks-products-and-mathworks-online-services-use)
