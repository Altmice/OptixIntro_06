# OptixIntro_06
## 1 基本公式
### 1.1 蒙特卡洛积分公式
$$\int_{a}^b{f(x)dx}=\frac{1}{N}\sum_{i=1}^N{\frac{f(X_i)}{p(X_i)}}\tag{1.1}$$
$$\int{\int{f(x,y)dxdy}}=\frac{1}{N}\sum_{i=1}^N{\frac{f(X_i,Y_i)}{p(X_i,Y_i)}}\tag{1.2}$$
### 1.2 采样分布公式
重要性采样指出，$p(X)$尽可能接近$f(X)$，而采样必须按照$p(X)$代表的分布进行。
#### 1.2.1 反演法
1、计算累积概率密度函数CDF $P(x)=\int_{0}^x{p(x')dx'}$。如果$p(X)$等于$f(X)$，反演法毫无意义。

2、计算反函数$P^{-1}(x)$

3、计算均匀分布随机变量$ξ$

4、计算$X_i=P^{-1}(ξ)$

相关资料参考*Siggraph03*的Arvo论文，*RaytracingGems*的第4章。

例如代码中的UnitSquareToCosineHemisphere，在半球空间是余弦分布$p(θ)$

$$\int_{Ω}{p(θ)dω}=\int_{0}^{φ}{\int_{0}^{θ}A\cos{θ'}\sin{θ'}dθ'dφ'}$$

将函数拆为两个积分

$$\int_{0}^{φ}{\frac{A}{2π}dφ'}=\frac{Aφ}{2π}=ξ_1$$
和

$$\int_{0}^{θ}2π\cos(θ')\sin(θ')dθ'=π(1-\cos^2(θ))=ξ_2$$

首先归一化概率，令$θ=π/2,φ=2π,\int_{Ω}{p(θ)dω}=1$,解得

$$A=\frac{1}{π}$$

取反函数$φ=2πξ_1,\cos(θ)=\sqrt{1-ξ_2}$,那么方向的三维分布是

$$(x,y,z)=(\sinθ\cosφ,\sinθ\sinφ,\cosθ)=(\cos(2πξ_1)\sqrt{ξ_2},\sin(2πξ_1)\sqrt{ξ_2},\sqrt{1-ξ_2})$$

随后可以按照正交轴方向变换过去

$$RandomVector=x*\boldsymbol{s}+y*\boldsymbol{t}+z*\boldsymbol{n}$$

#### 1.2.2 舍选法

假定$x\sim{p}$,$p:[a,b]\mapsto{\mathbb{R}}$,对于任意x，$p(x)<m$,由此我们在矩形$[a,b]\times{[0,m]}$采样，选择$y<p(x)$:
```c
done = false;
while(not done)
{
x=a+r()*(b-a);
y=r()*m;
if (y < p(x))
    done = true

}
```
#### 1.2.3   Metropolis
Metropolis的思想是从一个随机采样点出发，通过扰动变换到下一个采样点。

### 1.3 路径积分公式
Kajiya渲染公式，出射辐射度表示为物体表面发光和来自空间其他方向光源的入射的折反射之和。

$$L_o(p,\omega_o{})=L_e(p,\omega_o{})+\int_\mathbb{S^2}{f(p,\omega_o{},\omega_i{})L_d(p,\omega_i{})|\cos{θ_i}|d\omega_i}\tag{1.3}$$
### 1.4 多重重要性采样
可以看到，光照积分是光源和BSDF的乘积。对于$f(x)g(x)$的积分，当难以按照其乘积做重点采样，不妨分别按$p_f$和$p_g$采样。基于重要性采样的蒙特卡洛积分函数为

$$\frac{1}{n_f}\sum_{i=1}^{n_f}{\frac{f(X_i)g(X_i)w_f(X_i)}{p_f(X_i)}}+\frac{1}{n_g}\sum_{j=1}^{n_g}{\frac{f(X_j)g(X_j)w_g(X_j)}{p_g(X_j)}}\tag{1.4}$$
平衡自适应采样系数

$$w_s(s)=\frac{n_sp_s(x)}{\sum_i{n_ip_i(x)}}\tag{1.5}$$
幂自适应系数

$$w_s(x)=\frac{(n_sp_s(x))^\beta}{\sum_i{n_ip_i(x)}^\beta}\tag{1.6}$$
Veach指出$\beta=2$效果较好

OptixIntro_06给出了$n_f=n_g=1$的情形
```c
float balanceHeuristic(const float a,const float b)
{
    return a/(a+b);
}
float powerHeuristic(const float a,const float b)
{
    const float t=a*a;
    return t/(t+b*b);
}
```
### 1.5 轮盘赌和划分机制
按照Whitted光线分裂的追迹方法，适用于光学表面而不适用于粗糙面。
当函数是由一系列加和组成，比如光线分裂后的空间采样光线之和

$$F=F_1+…+F_N$$

考虑到相当一部分求和项的贡献很少，而跟踪的成本很高，所以加和项替换为下面的形式

$$F_i'=\lbrace\begin{matrix}\frac{F_i}{q_i}&ξ<q_i\\0&Otherwise\end{matrix}\tag{1.7}$$

注意$E[F_i']=\frac{E[F_i]}{q_i}*q_i=E[F_i]$，因此是一种无偏(unbiased)的方法。
划分机制指

$$F_i'=\frac{1}{k}\sum_{j=1}^k{F_{i,j}}\tag{1.8}$$
划分系数k根据估算贡献值选，贡献值高划分k越多。因此通过轮盘赌和划分，对于低贡献的路径尽早截止，对于高贡献的路径则丰富采样。

##  2 相机模型和积分器
光线追迹的开始，是依据特定的相机模型，生成光线。自定义的光线信息结构体为
```c
optix::float4 absorption_ior;//xyz是吸收系数，w是当前材料折射率
optix::float2 ior;//.x 光线在内部 .y 外折射率
optix::float3 pos;//当前世界坐标系下的落点
float distance;//吸收材料用的传播距离
optix::float3 wo,wi;//同一条光线，wo朝向观察者方向，wi朝向光源方向

optix::float3 radiance;//积分结果辐射度
int flags;//标记散射、镜面、体散射等等

optix::float3 f_over_pdf;//被积函数f(x)
float pdf;//采样概率分布p(x)

optix::float3 extinction;//体吸收系数
unsigned int seed;
```
OptixIntro_06提供了三种相机模型生成光线，分辨率为1spp。
随后对光线进行积分，并一帧帧地累加平均结果。
重点在于积分器函数
```c

RT_FUNCTION void integrator(PerRayData& prd, float3& radiance)
{
  // This renderer supports nested volumes. Four levels is plenty enough for most cases.渲染器支持夹心材质，安排了四层。
  // The absorption coefficient and IOR of the volume the ray is currently inside. RGB吸收系数，和光线所处材质的折射率
  float4 absorptionStack[MATERIAL_STACK_SIZE]; // .xyz == absorptionCoefficient (sigma_a), .w == index of refraction

  radiance = make_float3(0.0f); // Start with black.待计算的辐射度
  
  float3 throughput = make_float3(1.0f); // The throughput for the next radiance, starts with 1.0f. 下一级光线的系数

  int stackIdx = MATERIAL_STACK_EMPTY; // Start with empty nested materials stack.  -1 空气层
  int depth = 0;                       // Path segment index. Primary ray is 0. 用于限制追迹深度的关键变量

  prd.absorption_ior = make_float4(0.0f, 0.0f, 0.0f, 1.0f); // Assume primary ray starts in vacuum. 光线从真空开始，折射率1.0f，不存在吸收

  prd.flags = 0;

  // Russian Roulette path termination after a specified number of bounces needs the current depth.光线最大弹跳深度限制
  while (depth < sysPathLengths.y)
  {
    prd.wo        = -prd.wi;           // Direction to observer.向观察者(相机)方向
    prd.ior       = make_float2(1.0f); // Reset the volume IORs.重置折射率
    prd.distance  = RT_DEFAULT_MAX;    // Shoot the next ray with maximum length.
    prd.flags    &= FLAG_CLEAR_MASK;   // Clear all non-persistent flags. In this demo only the last diffuse surface interaction stays.清楚除了漫反射以外的所有标志

    // Handle volume absorption of nested materials.处理夹心材料吸收
    if (MATERIAL_STACK_FIRST <= stackIdx) // Inside a volume?
    {
      prd.flags     |= FLAG_VOLUME;                            // Indicate that we're inside a volume. => At least absorption calculation needs to happen.
      prd.extinction = make_float3(absorptionStack[stackIdx]); // There is only volume absorption in this demo, no volume scattering.取吸收系数
      prd.ior.x      = absorptionStack[stackIdx].w;            // The IOR of the volume we're inside. Needed for eta calculations in transparent materials.介质内折射率
      if (MATERIAL_STACK_FIRST <= stackIdx - 1)
      {
        prd.ior.y = absorptionStack[stackIdx - 1].w; // The IOR of the surrounding volume. Needed when potentially leaving a volume to calculate eta in transparent materials.外层的包围介质折射率
      }
    }

    // Note that the primary rays (or volume scattering miss cases) wouldn't noramlly offset the ray t_min by sysSceneEpsilon. Keep it simple here.
    optix::Ray ray = optix::make_Ray(prd.pos, prd.wi, 0, sysSceneEpsilon, prd.distance);
    rtTrace(sysTopObject, ray, prd);

    // This renderer supports nested volumes.夹层渲染
    if (prd.flags & FLAG_VOLUME)
    {
      // We're inside a volume. Calculate the extinction along the current path segment in any case.夹层渲染计算其沿路径的损耗吸收
      // The transmittance along the current path segment inside a volume needs to attenuate the ray throughput with the extinction
      // before it modulates the radiance of the hitpoint.
      throughput *= expf(-prd.distance * prd.extinction);
    }

    radiance += throughput * prd.radiance;//增加本层吸收后的辐射度

    // Path termination by miss shader or sample() routines.
    // If terminate is true, f_over_pdf and pdf might be undefined.光线已经丢失
    if ((prd.flags & FLAG_TERMINATE) || prd.pdf <= 0.0f || isNull(prd.f_over_pdf))
    {
      break;
    }

    // PERF f_over_pdf already contains the proper throughput adjustment for diffuse materials: f * (fabsf(optix::dot(prd.wi, state.normal)) / prd.pdf);  f(x)cos(θ)/p(x),throughput此时代表需要求和的Fi
    throughput *= prd.f_over_pdf;

    // Unbiased Russian Roulette path termination.
    if (sysPathLengths.x <= depth) // Start termination after a minimum number of bounces.几次弹跳后开始截止
    {
      const float probability = fmaxf(throughput); // DAR Other options: // intensity(throughput); // fminf(0.5f, intensity(throughput));
      if (probability < rng(prd.seed)) // Paths with lower probability to continue are terminated earlier.值太小了
      {
        break;
      }
      throughput /= probability; // Path isn't terminated. Adjust the throughput so that the average is right again.期望系数调整，按最大的归一化
    }

    // Adjust the material volume stack if the geometry is not thin-walled but a border between two volumes 处理夹心介质的进栈和出栈
    // and the outgoing ray direction was a transmission.
    if ((prd.flags & (FLAG_THINWALLED | FLAG_TRANSMISSION)) == FLAG_TRANSMISSION) 
    {
      // Transmission.
      if (prd.flags & FLAG_FRONTFACE) // Entered a new volume?
      {
        // Push the entered material's volume properties onto the volume stack.
        //rtAssert((stackIdx < MATERIAL_STACK_LAST), 1); // Overflow?
        stackIdx = min(stackIdx + 1, MATERIAL_STACK_LAST);
        absorptionStack[stackIdx] = prd.absorption_ior;
      }
      else // Exited the current volume?
      {
        // Pop the top of stack material volume.
        // This assert fires and is intended because I tuned the frontface checks so that there are more exits than enters at silhouettes.
        //rtAssert((MATERIAL_STACK_EMPTY < stackIdx), 0); // Underflow?
        stackIdx = max(stackIdx - 1, MATERIAL_STACK_EMPTY);
      }
    }

    ++depth; // Next path segment.
  }
}
```
## 3 光照计算
在[Detlef的视频](http://on-demand.gputechconf.com/gtc/2018/video/S8518/ "Introduction To Optix")中，提出了所谓的**USE_ NEXT_EVENT_ESTIMATION**，主要针对光线在粗糙面发生散射的情形，光线在粗糙面会主动地寻找光源计算亮度。一条长度为n的路径，其贡献值不是只有最终到光源或者丢失的部分，而是在每一级都尽可能寻找光源产生贡献值。这样图像就能以更快的速度收敛。

环境中的光照由三部分组成，追迹到光源的直接光照，在粗糙表面弹跳的间接光照和环境均匀光照。分别对应*closesthit_light.cu*, *closesthit.cu* 以及 *miss.cu*.我们重点分析前两者。

### 3.1 到达光源的直接光照
这一部分的主要问题是空间面光源对于散射表面的pdf，跟据[Ray Tracing: The rest of your life](https://github.com/petershirley/raytracingtherestofyourlife/blob/master/README.md "下载页面")第7章的推导方法

$$d\omega=dA\cos(\alpha)/Distance^2$$

$Distance$是光线传播距离，$\alpha$是光源表面法线和光线方向的夹角。
按照空间角均匀采样的概率和按面积采样的概率应该相等，所以

$$\frac{p(\omega)dA\cos(\alpha)}{L^2}=\frac{dA}{A}$$
$$p(\omega)=\frac{L^2}{A\cos{\alpha}}\tag{1.9}$$
一般来说，概率值会很大，因为相对的空间角很小。
```c
// Very simple closest hit program just for rectangle area lights.
// 用于矩形光源的Closehit程序
RT_PROGRAM void closesthit_light()
{
  thePrd.pos      = theRay.origin + theRay.direction * theIntersectionDistance; // Advance the path to the hit position in world coordinates. 光线落点
  thePrd.distance = theIntersectionDistance; // Return the current path segment distance, needed for absorption calculations in the integrator.传播距离，跟计算空间角有关

  const float3 geoNormal = optix::normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, varGeoNormal)); // PERF Not really needed when it's know that light geometry is not under Transforms.几何法线

  const float cosTheta = optix::dot(thePrd.wo, geoNormal);
  thePrd.flags |= (0.0f <= cosTheta) ? FLAG_FRONTFACE : 0;//确定前后方向

  thePrd.radiance = make_float3(0.0f); // Backside is black.

  if (thePrd.flags & FLAG_FRONTFACE) // Looking at the front face?
  {
    const LightDefinition light = sysLightDefinitions[parLightIndex];
    
    thePrd.radiance = light.emission;//取光源亮度

#if USE_NEXT_EVENT_ESTIMATION
    const float pdfLight = (thePrd.distance * thePrd.distance) / (light.area * cosTheta); // Solid angle pdf. Assumes light.area != 0.0f.注意是空间角的倒数，这个pdf基本很大
    // If it's an implicit light hit from a diffuse scattering event and the light emission was not returning a zero pdf.
    if ((thePrd.flags & FLAG_DIFFUSE) && DENOMINATOR_EPSILON < pdfLight)
    {
      //这里必须是散射面，散射面多重重要性采样的两部分，由于thePrd.pdf<<pdfLight，这里radiance基本是0.
      // Scale the emission with the power heuristic between the previous BSDF sample pdf and this implicit light sample pdf.
      thePrd.radiance *= powerHeuristic(thePrd.pdf, pdfLight);
    }
#endif // USE_NEXT_EVENT_ESTIMATION
  }

  // Lights have no other material properties than emission in this demo. Terminate the path.
  thePrd.flags |= FLAG_TERMINATE;
}
```
### 3.2 物体表面的间接光照

```c
RT_PROGRAM void closesthit()
{
  State state; // All in world space coordinates!

  state.geoNormal = optix::normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, varGeoNormal));
  state.normal    = optix::normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, varNormal));

  thePrd.pos      = theRay.origin + theRay.direction * theIntersectionDistance; // Advance the path to the hit position in world coordinates.
  thePrd.distance = theIntersectionDistance; // Return the current path segment distance, needed for absorption calculations in the integrator.

  // Explicitly include edge-on cases as frontface condition!
  // Keeps the material stack from overflowing at silhouttes.
  // Prevents that silhouettes of thin-walled materials use the backface material.
  // Using the true geometry normal attribute as originally defined on the frontface!
  thePrd.flags |= (0.0f <= optix::dot(thePrd.wo, state.geoNormal)) ? FLAG_FRONTFACE : 0;

  if ((thePrd.flags & FLAG_FRONTFACE) == 0) // Looking at the backface?
  {
    // Means geometric normal and shading normal are always defined on the side currently looked at.
    // This gives the backfaces of opaque BSDFs a defined result.
    state.geoNormal = -state.geoNormal;
    state.normal    = -state.normal;
    // Do not recalculate the frontface condition!
  }

  // A material system with support for arbitrary mesh lights would evaluate its emission here.
  // But since only parallelogram area lights are supported, those get a dedicated closest hit program to simplify this demo.
  thePrd.radiance = make_float3(0.0f);

  MaterialParameter parameters = sysMaterialParameters[parMaterialIndex];

  // Start fresh with the next BSDF sample.  (Either of these values remaining zero is an end-of-path condition.)
  thePrd.f_over_pdf = make_float3(0.0f);
  thePrd.pdf        = 0.0f;

  // Only the last diffuse hit is tracked for multiple importance sampling of implicit light hits.
  thePrd.flags = (thePrd.flags & ~FLAG_DIFFUSE) | parameters.flags; // FLAG_THINWALLED can be set directly from the material parameters.
//采样BSDF产生下一级方向wi,f_over_pdf(BSDF值)和pdf
  sysSampleBSDF[parameters.indexBSDF](parameters, state, thePrd);

#if USE_NEXT_EVENT_ESTIMATION
  // Direct lighting if the sampled BSDF was diffuse and any light is in the scene.
  if ((thePrd.flags & FLAG_DIFFUSE) && 0 < sysNumLights)
  {
    const float2 sample = rng2(thePrd.seed); // Use higher dimension samples for the position. (Irrelevant for the LCG).

    LightSample lightSample; // Sample one of many lights. 
  
    // The caller picks the light to sample. Make sure the index stays in the bounds of the sysLightDefinitions array.随机选一个光源
    lightSample.index = optix::clamp(static_cast<int>(floorf(rng(thePrd.seed) * sysNumLights)), 0, sysNumLights - 1); 

    const LightType lightType = sysLightDefinitions[lightSample.index].type;
//产生新的光线方向、概率、能量等等
    sysSampleLight[lightType](thePrd.pos, sample, lightSample); 
  
    if (0.0f < lightSample.pdf) // Useful light sample?
    {
      // Evaluate the BSDF in the light sample direction. Normally cheaper than shooting rays.
      // Returns BSDF f in .xyz and the BSDF pdf in .w
    //基于采样得到的方向，计算BSDF
      const float4 bsdf_pdf = sysEvalBSDF[parameters.indexBSDF](parameters, state, thePrd, lightSample.direction);

      if (0.0f < bsdf_pdf.w && isNotNull(make_float3(bsdf_pdf)))
      {
        // Do the visibility check of the light sample.
        PerRayData_shadow prdShadow;
      
        prdShadow.visible = true; // Initialize for miss.

        // Note that the sysSceneEpsilon is applied on both sides of the shadow ray [t_min, t_max] interval 
        // to prevent self intersections with the actual light geometry in the scene!
        optix::Ray ray = optix::make_Ray(thePrd.pos, lightSample.direction, 1, sysSceneEpsilon, lightSample.distance - sysSceneEpsilon); // Shadow ray.
        rtTrace(sysTopObject, ray, prdShadow);
//用阴影光线验证可见性
        if (prdShadow.visible)
        {
          if (thePrd.flags & FLAG_VOLUME) // Supporting nested materials includes having lights inside a volume.
          {
            // Calculate the transmittance along the light sample's distance in case it's inside a volume.
            // The light must be in the same volume or it would have been shadowed! 介质吸收
            lightSample.emission *= expf(-lightSample.distance * thePrd.extinction);
          }
//同样的，lightSample.pdf>>bsdf_pdf.w，misWeight接近于1
          const float misWeight = powerHeuristic(lightSample.pdf, bsdf_pdf.w);//多重重要性采样
          
          thePrd.radiance += make_float3(bsdf_pdf) * lightSample.emission * (misWeight * optix::dot(lightSample.direction, state.normal) / lightSample.pdf);
        }
      }
    }
  }
#endif // USE_NEXT_EVENT_ESTIMATION
}

```
 这里针对散射面的积分公式是

$$\frac{bsdf*L*w(pdf_L)*\cos{\alpha}}{pdf_L}$$
散射方向是跟据光源采样得到的，而不是bsdf。
多重重要性采样由两部分组成，按照Bsdf采样的另一部分放到closesthit_light()或者是miss_environment_constant()去了。对closesthit_light()来说，由于lightSample.pdf>>thePrd.pdf,因此这部分贡献接近于0。

## 4 方向采样
OptixIntro_06包含了三款材质，朗伯体(bsdf_diffuse_reflection.cu),镜面反射(bsdf_specular_reflection.cu)和(bsdf_specular_reflection_transmission)。
### 4.1 朗伯体材质
朗伯体$p(\omega)=\cos{\theta}/\pi$,方向按照反演法的示例计算。朗伯体表面是不产生能量的，因此

$$\frac{albedo*\cos{\theta}}{\pi*p(\omega)}=albedo$$

### 4.2 镜面反射材质
镜面反射是原有光路的延长。
```c
wi=reflect(-wo,normal);
f_over_pdf=albedo;
pdf=1.0f;
```
### 4.3 镜面透射材质
镜面透射是原光路的延长和分叉。这里主要用了Snell定律

$$n_i\sin{\theta_i}=n_t\sin{\theta_t}$$
和绝缘体菲涅尔折射率

$$r_\parallel=\frac{η_t\cos{\theta_i}-η_i\cos{\theta_t}}{η_t\cos{\theta_i}+η_i\cos{\theta_t}}$$

$$r_\perp=\frac{η_i\cos{\theta_i}-η_t\cos{\theta_t}}{η_i\cos{\theta_i}+η_t\cos{\theta_t}}$$
角标表示的是p偏振和s偏振，自然光反射率取两者的平均值，$t=1-r$,计算函数参考evaluateFresnelDielectric()

```c
计算反射方向和系数
计算透射方向和系数

按照反射系数做概率，随机选择反射或透射方向。此处可以理解成划分机制下的光线分裂，以及轮盘赌下的只计算一条光线代表全部。
相比Whitted，这样的处理避免了不断分叉带来的StackOverFlow。代价是需要花更多时间。

f_over_pdf=albedo;
pdf=1.0f;

```

4.2和4.3这里提一点，关于3.1提到的Distance，当光线从散射面经过光学表面落到光源时，Distance不应该仅仅是最后一面落点到光源点的距离。因为光学表面起到光路延长的作用，距离至少也应该从最后一个粗糙面点到光源点的传播距离算起。不过closesthit_light因为动用了幂次自适应系数，会导致在该函数中的贡献值很小。散射面的贡献获取主体在closesthit的对光源追迹。

本文没有把渲染方程推广到光线传输方程，看不懂的地方仍然存在。