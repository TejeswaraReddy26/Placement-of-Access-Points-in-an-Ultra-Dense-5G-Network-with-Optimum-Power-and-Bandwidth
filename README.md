 Placement of Access Points in an Ultra-Dense 5G Network with Optimum Power and Bandwidth

🚀 An intelligent system for optimizing the placement of 5G access points (APs) in ultra-dense networks using *Centroidal Voronoi Tessellation (CVT)* and *Machine Learning models* to achieve optimal coverage, reduced power consumption, and efficient bandwidth utilization.
## 📌 Overview
In ultra-dense 5G networks, traditional AP placement techniques (hexagonal/circular layouts) often lead to:
- ❌ Coverage overlap  
- ❌ Signal gaps  
- ❌ High power wastage  
This project proposes a *smart AP placement system* that:
- Ensures *uniform coverage*
- Reduces *power wastage (~15%)*
- Optimizes *bandwidth allocation*
- Supports *real-time prediction*
## 🧠 Core Idea

The system combines:

### 🔹 1. Centroidal Voronoi Tessellation (CVT)
- Places APs at *centroids of coverage regions*
- Minimizes overlap & uncovered zones
- Ensures balanced signal distribution

### 🔹 2. Machine Learning Models
- Logistic Regression (Best performer ~94% accuracy)
- SVM
- Random Forest
- Naive Bayes

- System Architecture

- <img width="1294" height="987" alt="image" src="https://github.com/user-attachments/assets/1d2799e6-c92a-46f4-8796-cb442d010216" />

Tech Stack

- 🐍 Python (Core ML + Backend)
- 📊 Scikit-learn (ML Models)
- 📈 NumPy / Pandas
- 🌐 Node.js (for real-time deployment)
- 📡 GIS / Signal Data Sources

## 🔄 Workflow

1. Collect raw data:
   - Geographic coordinates
   - Signal strength
   - User density

2. Preprocess data:
   - Remove noise
   - Normalize values

3. Extract features:
   - Signal strength
   - Location vectors
   - Demand patterns

4. Train ML models:
   - Logistic Regression (best accuracy)
   - Compare with SVM, RF, CVT

5. Predict optimal AP placement

6. Deploy for real-time use

## 📊 Results
| Model               | Accuracy |
|--------------------|---------|
| Logistic Regression | 94.2% ✅ |
| SVM                | 92.8% |
| CVT                | 93.1% |
| Random Forest      | 89.7% |
| Naive Bayes        | 87.5% |

📌 Logistic Regression provides the best balance between *accuracy & performance*.

## 🎯 Key Features
✔ Intelligent AP placement  
✔ Reduced power consumption (~15%)  
✔ Balanced coverage distribution  
✔ Real-time prediction capability  
✔ Scalable for smart cities  

## 🌍 Applications
- 🏙️ Smart Cities
- 📶 Urban 5G Deployment
- 🏢 Enterprise Networks
- 🚗 IoT & Connected Infrastructure

## ⚡ Performance
- ⚡ Predictions under *1 second*
- ⚡ Handles large datasets efficiently
- ⚡ Suitable for real-time deployment

## 🚧 Future Improvements
- 🔁 Reinforcement Learning for dynamic adaptation
- 🌦️ Real-time environmental data integration
- 🧠 Deep Learning models for complex scenarios
- 📊 Live dashboard visualization

## 📚 References
- IEEE WCNC (2018)
- 5G Network Research Papers
- Wireless Communication Journals
- 
## 👨‍💻 Authors
- Sandi Saketh Reddy  
- Nagireddy Subbareddy  
- Boilla Sai Tejeswara Reddy  
- Sandu Tharun  
- Dr. Pratham Majumder
- 
## ⭐ Contribute
Contributions are welcome!  
Feel free to fork, improve, and submit PRs 🚀

## 📜 License

This project is for academic and research purposes.
