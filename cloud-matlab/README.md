# Codex 云端 Linux MATLAB R2026a

当前目标：Debian 13、x86_64、2 核、8 GiB 内存，安装与用户本机一致的 MATLAB R2026a，使用已确认的学校账号许可。MATLAB 本体尚未安装；MathWorks 下载域名仍被云环境出站代理以 403 拒绝，账号登录与完整分析尚未实际验证。

安装路径默认是 `/workspace/MATLAB/R2026a`，可用 `MATLAB_ROOT` 指定其他路径。依赖清单采用 MathWorks 官方 R2026a/Debian 13 清单。缺失库及 `net-tools` 由 `runtime-deps.sh` 在无 root 权限下准备到 `/workspace/MATLAB/system-libs`，可用 `MATLAB_SYSTEM_LIBS` 修改位置。

## 安装前的环境条件

- 在云环境网络设置中放行 MathWorks 下载、Service Host 和账号认证需要的服务。至少包括 `www.mathworks.com`、`ssd.mathworks.com` 及其实际使用的下载/登录子域名；具体范围参考官方网络要求。脚本保留会话代理和 CA 校验，不能修改出站策略。
- 已确认账号关联 Campus-Wide 或 Individual 许可；实际可选的在线许可模式与运行权限仍需通过云端首次登录验证。本机登录状态不会自动迁移到云端。
- 账号模式需要适配 R2026a 的 MathWorks Service Host 及可写的用户运行数据/认证存储。当前用户目录只读；`MATLAB_PREFDIR` 只解决 MATLAB 偏好目录，不能视为已经解决 Service Host 的账号存储。需要平台提供受支持的可写用户数据环境，或由 MathWorks 确认其他受支持的配置。配置脚本自身仅管理工作目录，不改写 `HOME`，不搬运本机认证文件；MATLAB/Service Host 运行数据目录仍需平台解决。
- 当前安装盘约有 7.3 GiB、独立临时盘约有 4.9 GiB 可用。基础 MATLAB 约需 4.5 GB，安装峰值空间仍须实测；仅安装 MATLAB 本体并保留默认 JRE，临时文件放 `/tmp`。

## 安装与账号登录

以下命令在仓库根目录执行。默认版本为 R2026a，默认 `MATLAB_LICENSE_MODE=online`。

```bash
bash cloud-matlab/bootstrap.sh --check
bash cloud-matlab/bootstrap.sh --install-only
bash cloud-matlab/matlab.sh -nodisplay
```

第三条命令按官方 Login Named User 方式附加 `-licmode online`，供用户在交互终端完成首次登录。账号、密码和认证信息由用户在官方流程中处理，不写入代码或聊天。若许可证使用 New Licensing Technology，按 MathWorks 官方要求选择：

```bash
MATLAB_LICENSE_MODE=onlinelicensing bash cloud-matlab/matlab.sh -nodisplay
```

两个账号模式均忽略先前保存的网络许可端点，并为 MATLAB 子进程移除 `MLM_LICENSE_FILE`，避免选择错误的许可来源。该设置不会修改父进程环境或删除保存的端点。

账号登录与许可验证成功后，用同一模式运行：

```bash
bash cloud-matlab/matlab.sh -batch "disp(version)"
bash cloud-matlab/matlab.sh -batch "run('cloud-matlab/run_analysis.m')"
```

如果登录时选择了 `onlinelicensing`，在后续命令前同样加 `MATLAB_LICENSE_MODE=onlinelicensing`。首次登录与批处理复用是否成功均以真实 MATLAB 运行结果为准；目前尚未执行。

`--check` 只检查本地依赖和安装条件，不验证下载网络或许可签出。`--install-only` 安装文件后结束，不尝试登录。默认安装命令在账号模式下也不会自动进行首次账号认证。

## 网络许可证作为可选模式

只有明确选择网络许可时才使用以下方式；当前账号方案不要求取得学校 MLM 端口：

```bash
MLM_LICENSE_FILE='27000@license.example.edu' MATLAB_LICENSE_MODE=network bash cloud-matlab/matlab.sh -batch "disp(version)"
```

在本次会话，曾提供的网络许可端点保存于 `/workspace/work/matlab-bootstrap/license-server.txt`（权限 600）。仅当 `MATLAB_LICENSE_MODE=network` 且未显式设置 `MLM_LICENSE_FILE` 时读取；学校真实地址不提交到 GitHub。

网络模式需要服务器的 lmgrd 与 MLM 两个 TCP 端口都可达，可能需要学校 VPN。本平台尚未配置该 VPN/TCP 通道。普通 HTTP/HTTPS 代理不能代理 MATLAB 网络许可通信；仅写入端点不代表服务器可达。

## 运行及数据完整性

缺失依赖从 Debian 官方签名仓库下载，APT 索引写入工作目录，不修改系统 APT 状态。启动脚本为 MATLAB 添加用户目录中的库和工具路径，并将偏好目录设为 `/workspace/work/matlab-preferences/R2026a`；可用 `MATLAB_PREFDIR` 覆盖。

`run_analysis.m` 将 `genggai1101.m` 原样复制为 `/workspace/work/matlab-validation/fdm_paper_analysis.m`，使文件名与主函数一致。它生成图片但不显示窗口，输出至 `/workspace/work/matlab-validation/results`，不修改原统计方法。

运行成功后核查：全部组号 1–28；10CF 和 15CF 各 14 组；各材料的数据、实测排名和权重敏感性均包含相同的 14 个组号；排除组号为空；每项响应诊断 N=14；汇总仅两行 All 情景。

## 官方依据

- [R2026a Linux 系统需求](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2026a-linux.pdf)
- [R2026a / Debian 13 运行库清单](https://github.com/mathworks-ref-arch/container-images/blob/main/matlab-deps/r2026a/debian13/base-dependencies-amd64.txt)
- [MPM 安装参数、TMPDIR 和代理支持](https://www.mathworks.com/help/install/ug/mpminstall.html)
- [云端许可要求](https://www.mathworks.com/help/install/license/licensing-for-mathworks-products-running-on-the-cloud.html)
- [无显示器环境的账号登录方式](https://www.mathworks.com/matlabcentral/answers/547302-is-login-named-user-available-on-machines-without-a-display)
- [Login Named User 要求](https://www.mathworks.com/matlabcentral/answers/168441-what-is-login-named-user-and-how-can-i-use-it)
- [Service Host 管理和版本要求](https://github.com/mathworks-ref-arch/administer-mathworks-service-host)
- [MATLAB 偏好目录设置](https://www.mathworks.com/matlabcentral/answers/93696-how-do-i-change-the-matlab-preferences-directory-location)
- [MathWorks 网络要求](https://www.mathworks.com/matlabcentral/answers/97265-which-ports-and-addresses-do-mathworks-products-and-mathworks-online-services-use)
